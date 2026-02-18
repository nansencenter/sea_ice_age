from argparse import ArgumentParser
import datetime as dt
import glob
import os
from time import sleep
from string import Template
import subprocess
import hashlib

def run_subprocess(cmd, num_tries=1, sleep_time=20, stdout=subprocess.PIPE, shell=True, timeout=None):
    '''
    call subprocess.run with command and raise an error if it fails

    Parameters:
    -----------
    cmd : str
        command to be run through the system
    num_tries : int
        number of times to attempt the command
    sleep_time : int
        time to wait between retries (in seconds)
    shell : bool
        Create subprocess with shell emulation? Input to subprocess.run
    timeout : int
        Timeout of the subprocess
    stdout : int
        Output of the subprocess

    Returns:
    --------
    res : subprocess.CompletedProcess
        has attributes returncode, stdout, stderr
    '''
    for i in range(num_tries):
        # start subprocess in background
        
        res = subprocess.run(cmd,
                stdout=stdout,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                shell=shell,
                timeout=timeout,
                )
        if res.returncode == 0:
            # if success, return subprocess
            return res
        if i < num_tries - 1:
            #else, wait and try again
            print(res.stdout)
            print(f'Retrying in {sleep_time} seconds...')
            sleep(sleep_time)
    print(res.stdout)
    raise ValueError

product = "SEAICE_ARC_PHY_AUTO_L4_MY_011_025"
dataset = "cmems_obs-si_arc_phy-siage_my_P1D-m_202603"

xml_template = Template("""<?xml version='1.0' encoding='UTF-8'?><delivery product="$PRODUCT" PushingEntity="SIW-METNO-OSLO-NO" date="$SEND_DATE"> <dataset DatasetName="$DATASET">
$FILES
</dataset></delivery>""")

file_template = Template('<file FileName="$FILENAME" FileType="DT" FinalStatus="Delivered" Checksum="$CHECKSUM" StartUploadTime="$START_UPLOAD_TIME" StopUploadTime="$STOP_UPLOAD_TIME" />')


bookmark = 'cmems-my'
root_cmems_dir = 'NERSC_arctic25km_sea_ice_age_v2p1/zarr/cmems'


run_command = True
sleep_push_seconds = 0.1
sleep_response_seconds = 60
sleep_dir_seconds = 60

parser = ArgumentParser()
parser.add_argument('-im', "--input-dir-mask", type=str, help="Input directory mask", default="*/*")
parser.add_argument('-s', "--stop", action='store_true', help="Stop after the first directory (for testing)")

args = parser.parse_args()

input_dirs = sorted(glob.glob(f'{root_cmems_dir}/{args.input_dir_mask}'))
print(f"Found {len(input_dirs)} input directories:\n {input_dirs[0]}\n {input_dirs[-1]}")

for input_dir in input_dirs:
    date_dir = input_dir.replace(root_cmems_dir, '').lstrip('/')
    dnt_response_mask = f'{product}_{date_dir.replace("/", "_")}_*_response.xml'
    dnt_response_files = sorted(glob.glob(dnt_response_mask))
    dir_is_pushed = False
    if len(dnt_response_files) > 0:
        last_dnt_response_file = dnt_response_files[-1]
        with open(os.path.basename(last_dnt_response_file), 'r') as f:
            response_content = f.read()
        if 'Validated="True"' in response_content and 'Ingested="True"' in response_content:
            dir_is_pushed = True
    if dir_is_pushed:
        print(f'{input_dir} is already pushed, skipping')
        continue
    print(f'{input_dir} is not pushed, pushing now')

    netcdf_files = sorted(glob.glob(f'{input_dir}/*.nc'))
    print(f'Found {len(netcdf_files)} netCDF files in {input_dir} {netcdf_files[0]}')

    # first create the directory
    dst_dir = f'{product}/{dataset}/{date_dir}'
    lftp_cmd = f"mkdir -p {dst_dir}"
    cmd = f'lftp -e "{lftp_cmd}; bye;" {bookmark}'
    print(cmd)
    if run_command:
        run_subprocess(cmd, num_tries=10, timeout=1000)
        sleep(sleep_push_seconds)

    FILES = ''
    for filename in netcdf_files:
        # upload the file
        lftp_cmd = f"cd {dst_dir} && put {filename}"
        cmd = f'lftp -e "{lftp_cmd}; bye;" {bookmark}'
        t0 = dt.datetime.now(dt.UTC).replace(tzinfo=None)
        if run_command:
            run_subprocess(cmd, num_tries=10, timeout=1000)
            sleep(sleep_push_seconds)
        t1 = dt.datetime.now(dt.UTC).replace(tzinfo=None)

        # calculate the checksum
        with open(filename, "rb") as f:
            CHECKSUM = hashlib.md5(f.read()).hexdigest()

        # create the file XML content
        FILENAME = filename.replace(root_cmems_dir, '').lstrip('/')
        START_UPLOAD_TIME = t0.strftime('%Y%m%dT%H%M%S')
        STOP_UPLOAD_TIME = t1.strftime('%Y%m%dT%H%M%S')
        file_xml_content = file_template.substitute(FILENAME=FILENAME, CHECKSUM=CHECKSUM, START_UPLOAD_TIME=START_UPLOAD_TIME, STOP_UPLOAD_TIME=STOP_UPLOAD_TIME)
        FILES += file_xml_content + '\n'

    # create the DNT XML content
    SEND_DATE = dt.datetime.now(dt.UTC).strftime('%Y%m%dT%H%M%S')
    xml_content = xml_template.substitute(PRODUCT=product, DATASET=dataset, SEND_DATE=SEND_DATE, FILES=FILES)
    dnt_filename = f'{product}_{date_dir.replace("/", "_")}_{SEND_DATE}.xml'
    print(dnt_filename)
    with open(dnt_filename, 'w') as f:
        f.write(xml_content)

    # upload the DNT file
    dnt_dir = f'{product}/DNT'
    lftp_cmd = f"cd {dnt_dir} && put {dnt_filename}"
    cmd = f'lftp -e "{lftp_cmd}; bye;" {bookmark}'
    print(cmd)
    if run_command:
        run_subprocess(cmd, num_tries=10, timeout=1000)
        sleep(sleep_push_seconds)

    # check for the response file
    resp_dir = f'{product}/DNT_response'
    resp_filename = f'{resp_dir}/{dnt_filename.replace(".xml", "_response.xml")}'
    lftp_cmd = f"set xfer:clobber on; get {resp_filename}"
    cmd = f'lftp -e "{lftp_cmd}; bye;" {bookmark}'
    print(cmd)
    sleep(sleep_response_seconds)
    if run_command:
        run_subprocess(cmd, num_tries=10, timeout=1000)
        with open(os.path.basename(resp_filename), 'r') as f:
            response_content = f.read()
        if 'Validated="True"' not in response_content or 'Ingested="True"' not in response_content:
            raise ValueError("Response indicates validation failed")
        else:
            print(f"\n\n\t\tSUCESS: {input_dir}\n\n")

    if args.stop:
        break
    sleep(sleep_dir_seconds)
from argparse import ArgumentParser
import datetime as dt
import glob
import os
from time import sleep
from string import Template
import subprocess
import hashlib
from typing import List, Tuple


# Constants
PRODUCT = "SEAICE_ARC_PHY_AUTO_L4_MY_011_025"
DATASET = "cmems_obs-si_arc_phy-siage_my_P1D-m_202603"
BOOKMARK = 'cmems-my'
ROOT_CMEMS_DIR = 'NERSC_arctic25km_sea_ice_age_v2p1/zarr/cmems'

SLEEP_PUSH_SECONDS = 0.1
SLEEP_RESPONSE_SECONDS = 60
SLEEP_DIR_SECONDS = 60

NUM_RETRIES = 10
TIMEOUT = 1000

# XML Templates
XML_TEMPLATE = Template(
    """<?xml version='1.0' encoding='UTF-8'?>"""
    """<delivery product="$PRODUCT" PushingEntity="SIW-METNO-OSLO-NO" date="$SEND_DATE">"""
    """ <dataset DatasetName="$DATASET">
$FILES
</dataset></delivery>"""
)

FILE_TEMPLATE = Template(
    '<file FileName="$FILENAME" FileType="DT" FinalStatus="Delivered" '
    'Checksum="$CHECKSUM" StartUploadTime="$START_UPLOAD_TIME" '
    'StopUploadTime="$STOP_UPLOAD_TIME" />'
)


class SubprocessError(Exception):
    """Exception raised when subprocess execution fails."""
    pass


def run_subprocess(cmd: str, num_tries: int = 1, sleep_time: int = 20, 
                   stdout=subprocess.PIPE, shell: bool = True, timeout: int = None):
    """
    Execute a subprocess command with retry logic.

    Parameters:
    -----------
    cmd : str
        Command to be run through the system
    num_tries : int
        Number of times to attempt the command
    sleep_time : int
        Time to wait between retries (in seconds)
    shell : bool
        Create subprocess with shell emulation
    timeout : int
        Timeout of the subprocess in seconds
    stdout : int
        Output of the subprocess

    Returns:
    --------
    res : subprocess.CompletedProcess
        Has attributes returncode, stdout, stderr

    Raises:
    -------
    SubprocessError
        If command fails after all retries
    """
    for i in range(num_tries):
        res = subprocess.run(
            cmd,
            stdout=stdout,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=shell,
            timeout=timeout,
        )
        if res.returncode == 0:
            return res
        
        if i < num_tries - 1:
            print(res.stdout)
            print(f'Retrying in {sleep_time} seconds...')
            sleep(sleep_time)
    
    print(res.stdout)
    raise SubprocessError(f"Command failed after {num_tries} attempts: {cmd}")


def run_lftp_command(lftp_cmd: str, bookmark: str = BOOKMARK, 
                     num_tries: int = NUM_RETRIES, timeout: int = TIMEOUT):
    """Execute an lftp command."""
    cmd = f'lftp -e "{lftp_cmd}; bye;" {bookmark}'
    print(f"Executing: {cmd}")
    run_subprocess(cmd, num_tries=num_tries, timeout=timeout)


def check_directory_pushed(date_dir: str, product: str = PRODUCT) -> bool:
    """
    Check if a directory has already been successfully pushed.
    
    Parameters:
    -----------
    date_dir : str
        The date directory path
    product : str
        The product name
        
    Returns:
    --------
    bool
        True if directory is already pushed and validated
    """
    dnt_response_mask = f'DNT_response/{product}_{date_dir.replace("/", "_")}_*_response.xml'
    dnt_response_files = sorted(glob.glob(dnt_response_mask))
    
    if len(dnt_response_files) == 0:
        return False
    
    last_dnt_response_file = dnt_response_files[-1]
    try:
        with open(last_dnt_response_file, 'r') as f:
            response_content = f.read()
        return 'Validated="True"' in response_content and 'Ingested="True"' in response_content
    except FileNotFoundError:
        return False


def calculate_checksum(filepath: str) -> str:
    """Calculate MD5 checksum of a file."""
    with open(filepath, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def create_remote_directory(dst_dir: str):
    """Create a remote directory using lftp."""
    lftp_cmd = f"cd {dst_dir} || mkdir -p {dst_dir}"
    run_lftp_command(lftp_cmd)
    sleep(SLEEP_PUSH_SECONDS)


def upload_file(filepath: str, dst_dir: str) -> Tuple[dt.datetime, dt.datetime]:
    """
    Upload a file to remote directory.
    
    Returns:
    --------
    Tuple[dt.datetime, dt.datetime]
        Start and stop upload times
    """
    lftp_cmd = f"cd {dst_dir} && put {filepath}"
    t0 = dt.datetime.now(dt.UTC).replace(tzinfo=None)
    run_lftp_command(lftp_cmd)
    sleep(SLEEP_PUSH_SECONDS)
    t1 = dt.datetime.now(dt.UTC).replace(tzinfo=None)
    return t0, t1


def create_file_xml_entry(filename: str, checksum: str, 
                          start_time: dt.datetime, stop_time: dt.datetime) -> str:
    """Create XML entry for a single file."""
    return FILE_TEMPLATE.substitute(
        FILENAME=filename,
        CHECKSUM=checksum,
        START_UPLOAD_TIME=start_time.strftime('%Y%m%dT%H%M%S'),
        STOP_UPLOAD_TIME=stop_time.strftime('%Y%m%dT%H%M%S')
    )


def create_dnt_xml(files_xml: str, product: str = PRODUCT, dataset: str = DATASET) -> str:
    """Create DNT XML content."""
    send_date = dt.datetime.now(dt.UTC).strftime('%Y%m%dT%H%M%S')
    return XML_TEMPLATE.substitute(
        PRODUCT=product,
        DATASET=dataset,
        SEND_DATE=send_date,
        FILES=files_xml
    )


def upload_dnt_file(dnt_filename: str, product: str = PRODUCT):
    """Upload DNT file to remote server."""
    remote_dnt_dir = f'{product}/DNT'
    local_dnt_path = f'DNT/{dnt_filename}'
    lftp_cmd = f"cd {remote_dnt_dir} && put {local_dnt_path}"
    run_lftp_command(lftp_cmd)
    sleep(SLEEP_PUSH_SECONDS)


def validate_response(resp_filename: str) -> bool:
    """
    Check and validate response file.
    
    Returns:
    --------
    bool
        True if validation and ingestion were successful
        
    Raises:
    -------
    ValueError
        If response indicates validation failed
    """
    local_resp_path = f'DNT_response/{os.path.basename(resp_filename)}'
    with open(local_resp_path, 'r') as f:
        response_content = f.read()
    
    validated = 'Validated="True"' in response_content
    ingested = 'Ingested="True"' in response_content
    
    if not validated or not ingested:
        raise ValueError(f"Response indicates validation failed: {response_content}")
    
    return True


def download_response_file(dnt_filename: str, product: str = PRODUCT):
    """Download and validate response file."""
    remote_resp_dir = f'{product}/DNT_response'
    remote_resp_filename = f'{remote_resp_dir}/{dnt_filename.replace(".xml", "_response.xml")}'
    local_resp_filename = os.path.basename(remote_resp_filename)
    lftp_cmd = f"set xfer:clobber on; cd {remote_resp_dir} && get {local_resp_filename} -o DNT_response/{local_resp_filename}"
    
    print(f"Waiting {SLEEP_RESPONSE_SECONDS}s for response file...")
    sleep(SLEEP_RESPONSE_SECONDS)
    
    run_lftp_command(lftp_cmd)
    validate_response(remote_resp_filename)


def process_netcdf_files(netcdf_files: List[str], dst_dir: str) -> str:
    """
    Process and upload NetCDF files, returning XML content.
    
    Parameters:
    -----------
    netcdf_files : List[str]
        List of NetCDF file paths to process
    dst_dir : str
        Destination directory on remote server
        
    Returns:
    --------
    str
        XML content for all uploaded files
    """
    files_xml = []
    
    for filepath in netcdf_files:
        # Upload file
        t0, t1 = upload_file(filepath, dst_dir)
        
        # Calculate checksum
        checksum = calculate_checksum(filepath)
        
        # Create XML entry
        filename = filepath.replace(ROOT_CMEMS_DIR, '').lstrip('/')
        file_xml = create_file_xml_entry(filename, checksum, t0, t1)
        files_xml.append(file_xml)
    
    return '\n'.join(files_xml) + '\n'


def process_directory(input_dir: str) -> bool:
    """
    Process a single directory: upload files and create DNT.
    
    Parameters:
    -----------
    input_dir : str
        Input directory path containing NetCDF files
        
    Returns:
    --------
    bool
        True if directory was actually processed, False if skipped
    """
    date_dir = input_dir.replace(ROOT_CMEMS_DIR, '').lstrip('/')
    
    # Check if already pushed
    if check_directory_pushed(date_dir):
        print(f'{input_dir} is already pushed, skipping')
        return False
    
    print(f'{input_dir} is not pushed, pushing now')
    
    # Find NetCDF files
    netcdf_files = sorted(glob.glob(f'{input_dir}/*.nc'))
    if not netcdf_files:
        print(f'No NetCDF files found in {input_dir}, skipping')
        return False
    
    print(f'Found {len(netcdf_files)} netCDF files in {input_dir}')
    print(f'  First: {netcdf_files[0]}')
    
    # Create local subdirectories if they don't exist
    os.makedirs('DNT', exist_ok=True)
    os.makedirs('DNT_response', exist_ok=True)
    
    # Create remote directory
    dst_dir = f'{PRODUCT}/{DATASET}/{date_dir}'
    create_remote_directory(dst_dir)
    
    # Process and upload files
    files_xml = process_netcdf_files(netcdf_files, dst_dir)
    
    # Create and save DNT XML
    xml_content = create_dnt_xml(files_xml)
    send_date = dt.datetime.now(dt.UTC).strftime('%Y%m%dT%H%M%S')
    dnt_filename = f'{PRODUCT}_{date_dir.replace("/", "_")}_{send_date}.xml'
    
    print(f'Creating DNT file: {dnt_filename}')
    with open(f'DNT/{dnt_filename}', 'w') as f:
        f.write(xml_content)
    
    # Upload DNT file
    upload_dnt_file(dnt_filename)
    
    # Download and validate response
    download_response_file(dnt_filename)
    
    print(f"\n\n\t\tSUCCESS: {input_dir}\n\n")
    return True


def main():
    """Main execution function."""
    parser = ArgumentParser(description="Push sea ice age data to CMEMS")
    parser.add_argument(
        '-im', "--input-dir-mask", 
        type=str, 
        help="Input directory mask (e.g., '2026/03' or '*/*')", 
        default="*/*"
    )
    parser.add_argument(
        '-s', "--stop", 
        action='store_true', 
        help="Stop after the first directory (for testing)"
    )

    args = parser.parse_args()
    
    # Find input directories
    input_dirs = sorted(glob.glob(f'{ROOT_CMEMS_DIR}/{args.input_dir_mask}'))
    
    if not input_dirs:
        print(f"No directories found matching: {ROOT_CMEMS_DIR}/{args.input_dir_mask}")
        return
    
    print(f"Found {len(input_dirs)} input directories:")
    print(f"  First: {input_dirs[0]}")
    print(f"  Last:  {input_dirs[-1]}")
    print()
    
    # Process each directory
    for i, input_dir in enumerate(input_dirs, 1):
        print(f"[{i}/{len(input_dirs)}] Processing: {input_dir}")
        
        try:
            was_processed = process_directory(input_dir)
        except Exception as e:
            print(f"ERROR processing {input_dir}: {e}")
            if args.stop:
                raise
            continue
        
        if args.stop and was_processed:
            print(f"Stopping after processing {input_dir} (--stop flag)")
            break
        
        if i < len(input_dirs) and was_processed:
            print(f"Waiting {SLEEP_DIR_SECONDS}s before next directory...")
            sleep(SLEEP_DIR_SECONDS)
    
    print(f"\nCompleted processing {len(input_dirs)} directories")


if __name__ == "__main__":
    main()
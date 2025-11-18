import os

import pandas as pd
import zarr
from zarr.storage import ZipStore

from .zarr_index import init_db, get_session, File, Array

class MeshFileHdf5:
    def __init__(self, filepath):
        self.filepath = filepath

    def read_names(self):
        if not os.path.exists(self.filepath):
            return []
        with pd.HDFStore(self.filepath, mode='r') as store:
            return list(store.keys())

    def load(self, read_names=()):
        data = {}
        if not read_names:
            read_names = self.read_names()
        for name in read_names:
            data[name] = pd.read_hdf(self.filepath, key=name)
        return data

    def save(self, data, mode='w'):
        for name in data:
            data[name].to_hdf(self.filepath, key=name, mode=mode)
            mode = 'a'


class MeshFileZar:
    """ Major problem with all varibale in one zip file is that
    if one variable is corrupted the whole file becomes unreadable.
    To fix corrupted zip file one needs to remove corrupted variable
    from the zip archive. But this is not straightforward with zarr API.
    """
    def __init__(self, filepath):
        self.filepath = filepath

    def read_names(self):
        if not os.path.exists(self.filepath):
            return []
        with ZipStore(self.filepath, mode='r') as store:
            z = zarr.open(store, mode='r')
            return list(z.array_keys())

    def load(self, read_names=(), as_dict=True):
        data = {}
        with ZipStore(self.filepath, mode='r') as store:
            z = zarr.open(store, mode='r')
            if read_names:
                z_array_keys = read_names
            else:
                z_array_keys = z.array_keys()
            for name in z_array_keys:
                data[name] = z[name][:]
        if not as_dict and read_names:
            data = [data[name] for name in read_names]
        return data

    def save(self, data, mode='w'):
        """ Save data to a Zarr file in ZIP format 
        Parameters
        ----------
        data : dict
            Dictionary where keys are names of arrays and values are numpy arrays to be saved.
        mode : str, optional, default = 'w
            'w': Write only the items in data to the file. Existing data in file is removed.
            'a': Append the items in data to the file. Existing data in file is preserved.
            'o': Overwrite existing items in the file with those in data. Existing items not in data are preserved.
        """
        if len(self.read_names()) == 0:
            mode = 'w'
        if mode == 'o':
            raw_data = self.load()
            raw_data.update(data)
            data = raw_data
            mode = 'w'
        with ZipStore(self.filepath, mode=mode) as store:
            z = zarr.group(store=store)
            for name, item in data.items():
                if name in z:
                    continue
                xa = z.create_array(name, shape=item.shape, dtype=item.dtype)
                xa[:] = item
        self.update_db()

    def remove(self, names):
        """ Remove arrays with specified names from the Zarr file.
        Parameters
        ----------
        names : list
            List of names of arrays to be removed.
        """
        data = self.load()
        data = {i: data[i] for i in data if i not in names}
        self.save(data, mode='w')

    def update_db(self):
        """ Update database index for this file """
        init_db()
        abs_path = os.path.abspath(self.filepath)
        st = os.stat(abs_path)
        mtime = st.st_mtime

        # current array names in the file (after write)
        array_names = self.read_names()
        with get_session() as session:
            db_file = session.get(File, abs_path)
            if db_file is None:
                db_file = File(path=abs_path, mtime=mtime)
                session.add(db_file)
            else:
                db_file.mtime = mtime

            # replace arrays entries for this file
            session.query(Array).filter(Array.path == abs_path).delete(synchronize_session=False)
            session.flush() # ensure deletion is executed before adding new entries
            for name in array_names:
                db_file.arrays.append(Array(path=abs_path, array_name=name))


class MeshFile(MeshFileZar):
    pass
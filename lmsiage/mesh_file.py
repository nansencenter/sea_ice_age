import os

import pandas as pd
import zarr
from zarr.storage import ZipStore

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
                

class MeshFile(MeshFileZar):
    pass
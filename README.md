# Installation
```
conda create -y -n sia gdal cartopy pip
conda activate sia
pip install pythesint netcdf4 nansat pygrib
python -c 'import pythesint as pti; pti.update_all_vocabularies()'
```


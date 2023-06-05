---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Dandiset Manifest
### This notebook shows all of the assets within this dandiset, organized by their transgenic lines and recording locations

```{code-cell} ipython3
from dandi.dandiapi import DandiAPIClient
```

```{code-cell} ipython3
import numpy as np
import pandas as pd
import pynwb
from tqdm import tqdm
from joblib import Parallel, delayed

import warnings
```

```{code-cell} ipython3
dandiset_id = '000039'
```

```{code-cell} ipython3
client = DandiAPIClient()
dandisets = list(client.get_dandisets())
```

```{code-cell} ipython3
ds = client.get_dandiset(dandiset_id)
```

```{code-cell} ipython3
asset_list = list(ds.get_assets())
number_assets = len(asset_list)
manifest = pd.DataFrame(columns=('session_id','specimen_id','genotype','area','imaging_depth','sex','age','path','date'), 
                        index=range(number_assets))

def to_struct(asset):
    manifest = {}
    md = asset.get_raw_metadata()
    manifest['session_id'] = md['wasGeneratedBy'][0]['identifier']
    manifest['specimen_id'] = md['wasAttributedTo'][0]['identifier']
    manifest['genotype'] = md['wasAttributedTo'][0]['genotype']
    manifest['sex'] = md['wasAttributedTo'][0]['sex']['name']
    manifest['age'] = md['wasAttributedTo'][0]['age']['value']
    manifest['path'] = md['path']
    manifest['date'] = md['wasGeneratedBy'][0]['startDate']
    
    path = md['path']
    s3_url = asset.get_content_url(regex='s3')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        io = pynwb.NWBHDF5IO(s3_url, mode='r', load_namespaces=True, driver='ros3')
        nwbfile = io.read()
    location = nwbfile.imaging_planes['imaging_plane_1'].location
    manifest['area'] = location.split(',')[0].split((' '))[1]
    manifest['imaging_depth'] = location.split(',')[1].split((' '))[1]
    return manifest

result = Parallel(n_jobs=10)(delayed(to_struct)(asset_list[i]) for i in tqdm(range(number_assets)))
```

```{code-cell} ipython3
manifest = pd.DataFrame(result)
```

```{code-cell} ipython3
manifest.head()
```

This dataframe has a row for each asset, describing some key metadata about the animal and recording location. We can explore the dataframe to identify assets by genotype, recording location (i.e. area), or other parameters, and find the path for each asset.

```{code-cell} ipython3
manifest.genotype.unique()
```

```{code-cell} ipython3
manifest.area.unique()
```

```{code-cell} ipython3
len(manifest[manifest.sex=='Female'])
```

```{code-cell} ipython3
pd.pivot_table(manifest, values=['path'],columns=['area'], index=['genotype'], aggfunc='count', fill_value=0)
```

```{code-cell} ipython3

```

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:26:15 2022

@author: lzhang481
"""

from dandi.dandiapi import DandiAPIClient
from pynwb import NWBHDF5IO
import pynwb
dandiset_id = '000021'  # ephys dataset from the Svoboda Lab
#filepath = 'sub-anm372795/sub-anm372795_ses-20170718.nwb'  # 450 kB file


with DandiAPIClient() as client:
asset = client.get_dandiset(dandiset_id, 'draft')
#get_asset_by_path(filepath)
#    s3_url = asset.get_content_url(follow_redirects=1, strip_query=True)
    

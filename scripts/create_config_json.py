#!/home/mraj/anaconda3/envs/Rv4/bin/python

"""
Script to be called from R to generate and save outer config file.

Started by Mukund Raj on 2022-01-25
"""

import json
import os
import re
import sys

arggs = str(sys.argv);
arggs = arggs.split(',')


# create info fields for outer config file from input arg string
zarr_store_name = re.sub(r'[^A-Za-z_.]+', '', arggs[1])
analysis_name = os.path.splitext(zarr_store_name)[0]
qc_object_name = analysis_name+'.qc'


# populate vitessce config

## read vitessce json template and customize the data

with open('../templates/view-config-inph-zarr.json') as json_file:
    data = json.load(json_file)


data["name"] = analysis_name
data_file_full_path_initial =  data["datasets"][0]["files"][0]["url"]
data_file_full_path = data_file_full_path_initial.replace("my_store", analysis_name)

data["datasets"][0]["files"][0]["url"] = data_file_full_path
data["datasets"][0]["files"][1]["url"] = data_file_full_path
data["datasets"][0]["files"][2]["url"] = data_file_full_path


vitessce_config = data

# populate outer config
outer_config = {}
outer_config["analysis_name"] = analysis_name
outer_config["zarr_store_name"] = zarr_store_name
outer_config["qc_object_name"] = qc_object_name
outer_config["vitessce_config"] = vitessce_config

# print outer_config

json_string = json.dumps(outer_config)

# Write using a JSON string
with open('../output/'+analysis_name+'.json', 'w') as outfile:
    outfile.write(json_string)

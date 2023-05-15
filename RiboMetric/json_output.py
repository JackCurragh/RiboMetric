"""
This script generates a machine readable format of the output in a json file.
This output can later be used to generate reports without needing to run the
parsing and modules again. 
"""
# Maybe json and csv in one script? Add option for keeping or leaving config? 
# Option to then use config in json (if exists) or local config
import json


def generate_json(
    results_dict: dict,
    config: dict,
    name: str = "RibosomeProfiler_data.json",
    outdir: str = "",
):
    """
    Generate a machine readable format of the RibosomeProfiler results

    Input:
        results_dict:
        config:
        name:
        outdir:

    Output:
        Writes to a json file
    """
    if "sequence_slice" in results_dict:
        del results_dict["sequence_slice"]
    
    if outdir == "":
        output = name
    else:
        if outdir.endswith("/") and outdir != "":
            outdir = outdir[:-1]
        output = outdir + "/" + name

    data = {"results": results_dict, "config": config}

    with open(output, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Results written in {output}")
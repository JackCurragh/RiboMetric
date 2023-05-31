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
    name: str = "RiboMetric_data.json",
    output_directory: str = "",
):
    """
    Generate a machine readable format of the RiboMetric results

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes to a json file
    """
    if "sequence_slice" in results_dict:
        del results_dict["sequence_slice"]

    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name + ".json"

    data = {"results": results_dict, "config": config}

    with open(output, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Results written in {output}")

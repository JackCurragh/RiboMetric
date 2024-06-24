"""
This script generates a machine readable format of the output in a json file.
This output can later be used to generate reports without needing to run the
parsing and modules again.
"""
# Maybe json and csv in one script? Add option for keeping or leaving config?
# Option to then use config in json (if exists) or local config
import json
import csv


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


def normalise_score(score: float, min_score: float, max_score: float) -> float:
    """
    Normalise the score of a metric

    Input:
        score: The score of the metric
        min_score: The minimum score of the metric
        max_score: The maximum score of the metric


    Output:
        The normalised score
    """
    if score > max_score:
        print(f"Score {score} is greater than the maximum score {max_score}")
        return 1
    elif score < min_score:
        return 0
    else:
        return (score - min_score) / (max_score - min_score)


def generate_csv(
    results_dict: dict,
    config: dict,
    name: str = "RiboMetric_data.json",
    output_directory: str = "",
):
    """
    Generate a csv file containing the different metrics and their
    corresponding score

    Input:
        results_dict: Dictionary containing the results of the qc analysis
        config: Dictionary containing the configuration information
        name: Name of the output file
        output_directory: Directory to write the output file to

    Output:
        Writes to a csv file
    """
    if output_directory == "":
        output = name
    else:
        if output_directory.endswith("/") and output_directory != "":
            output_directory = output_directory[:-1]
        output = output_directory + "/" + name + ".csv"

    columns = ["Metric", "Score", "MaxMinScore"]
    metrics_dict = []
    for key, value in results_dict["metrics"].items():
        if isinstance(value, float) or isinstance(value, int):
            if key not in config["max_mins"]:
                max_min_score = value
            else:
                max_min_score = normalise_score(
                    value,
                    config["max_mins"]['_'.join(key.split("_")[:-1])][0],
                    config["max_mins"]['_'.join(key.split("_")[:-1])][1]
                )
            metrics_dict.append(
                {"Metric": key,
                 "Score": value,
                 "MaxMinScore": max_min_score,
                 }
                 )
        elif isinstance(value, dict):
            for k, v in value.items():
                if key not in config["max_mins"]:
                    max_min_score = v
                else:
                    max_min_score  = normalise_score(
                        v,
                        config["max_mins"][key][0],
                        config["max_mins"][key][1]
                    )
                metrics_dict.append(
                    {
                        "Metric": f"{key}_{k}",
                        "Score": v,
                        "MaxMinScore": max_min_score
                        }
                        )

    with open(output, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        for data in metrics_dict:
            writer.writerow(data)
    print(f"Metrics written in {output}")

"""
Code in this script is used to generate the HTML and pdf output report
The functions are called by the main script RiboMetric.py
if the user specifies the --html flag
"""

from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from .modules import convert_html_to_pdf
import base64
import json

import os


def generate_report(
    plots: list,
    config: dict,
    export_mode: str = "html",
    name: str = "RiboMetric_report",
    outdir: str = "",
):
    """
    Generates a report of the RiboMetric results with plots

    Inputs:
        plots: A list containing the plots and metrics for the report
        export_mode: A string defining the mode of export: 'html', 'pdf' or
        'both' (Default: 'html')
        name: A string for the file name (Default: 'RiboMetric_report')
        outdir: A string for the output directory (Default: '')

    Outputs:
        No variables will be output
    """
    project_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = Environment(
        loader=FileSystemLoader(
            ["templates", f"{project_path}/RiboMetric/templates"]
        ),
        autoescape=False,
    )

    file_names = {"bam": config["argument"]["bam"].split("/")[-1]}
    if config["argument"]["annotation"] is not None:
        file_names["annotation"] = (config["argument"]["annotation"]
                                    .split("/")[-1])
    elif config["argument"]["gff"] is not None:
        file_names["annotation"] = config["argument"]["gff"].split("/")[-1]

    completion_time = datetime.now().strftime("%H:%M:%S %d/%m/%Y")

    binary_logo = open(
        f"{project_path}/RiboMetric/templates/RiboMetric_logo.png",
        "rb"
        ).read()
    base64_logo = base64.b64encode(binary_logo).decode("utf-8")

    binary_icon = open(
        f"{project_path}/RiboMetric/templates/RiboMetric_favicon.png",
        "rb"
        ).read()
    base64_icon = base64.b64encode(binary_icon).decode("utf-8")

    if outdir == "":
        output = name
    else:
        if outdir.endswith("/") and outdir != "":
            outdir = outdir[:-1]
        output = outdir + "/" + name

    if export_mode == "both":
        export_mode = ["html", "pdf"]
    else:
        export_mode = [export_mode]

    template = env.get_template("base.html")
    context = {
        "summary": plots.pop(0),
        "plots": plots,
        "file_names": file_names,
        "completion_time": completion_time,
        "logo": base64_logo,
        "favicon": base64_icon,
    }

    for filetype in export_mode:
        if filetype == "html":
            context["filetype"] = filetype
            jinja_render = template.render(context)
            out = output + ".html"
            with open(out, mode="w", encoding="utf-8") as f:
                f.write(jinja_render)
            print(f"Your {filetype} report can be found in {out}")
        else:
            context["filetype"] = filetype
            jinja_render = template.render(context)
            out = output + ".pdf"
            convert_html_to_pdf(jinja_render, out)
            print(f"Your {filetype} report can be found in {out}")


def int_keys_hook(data):
    """
    Custom object_hook for JSON parsing that converts number strings into
    integers
    """
    for key in list(data.keys()):
        if isinstance(key, str) and key.isdigit():
            data[int(key)] = data.pop(key)
    return data


def parse_json_input(json_path: str) -> tuple:
    """
    Parse json input from a previous RiboMetric run for use in generating plots

    Inputs:
        json_path: Path to json file

    Outputs:
        results_dict: Dictionary containing results from a RiboMetric analysis
        json_config: Config from the RiboMetric analysis
    """
    with open(json_path, 'r') as json_file:
        json_dict = json.load(json_file, object_hook=int_keys_hook)
    result_dict = json_dict["results"]
    json_config = json_dict["config"]
    return (result_dict, json_config)

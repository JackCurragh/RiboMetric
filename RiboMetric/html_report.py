"""
Code in this script is used to generate the HTML and pdf output report
The functions are called by the main script RiboMetric.py
if the user specifies the --html flag
"""

from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from .modules import convert_html_to_pdf
import base64

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
        results: A dictionary containing the results of
                    the RiboMetric analysis
        export_mode: A string defining the mode of
                    export: 'html','pdf' or 'both' (Default: 'html')
        name: A string for the file name (Default: 'RiboMetric_report')
        outdir: A string for the output directory (Default: '')

    Outputs:
        No variables will be output
    """
    env = Environment(
        loader=FileSystemLoader(["templates", "RiboMetric/templates"]),
        autoescape=False,
    )

    completion_time = datetime.now().strftime("%H:%M:%S %d/%m/%Y")

    binary_logo = open("RiboMetric/templates/RiboMetric_logo.png",
                       "rb").read()
    base64_logo = base64.b64encode(binary_logo).decode("utf-8")
    binary_icon = open("RiboMetric/templates/RiboMetric_favicon.png",
                       "rb").read()
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
        else:
            context["filetype"] = filetype
            jinja_render = template.render(context)
            out = output + ".pdf"
            convert_html_to_pdf(jinja_render, out)

    print(f"Done! Your report can be found in {out}")

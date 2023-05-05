"""
Code in this script is used to generate the HTML and pdf output report
The functions are called by the main script RibosomeProfiler.py
if the user specifies the --html flag
"""

from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from .modules import convert_html_to_pdf


def generate_report(
    plots: dict,
    export_mode: str = "html",
    name: str = "RibosomeProfiler_report",
    outdir: str = "",
):
    """
    Generates a report of the RibosomeProfiler results with plots

    Inputs:
        results: A dictionary containing the results of
                    the RibosomeProfiler analysis
        export_mode: A string defining the mode of
                    export: 'html','pdf' or 'both' (Default: 'html')
        name: A string for the file name (Default: 'RibosomeProfiler_report')
        outdir: A string for the output directory (Default: '')

    Outputs:
        No variables will be output
    """
    env = Environment(loader=FileSystemLoader("templates"), autoescape=False)

    completion_time = datetime.now().strftime("%H:%M:%S %d/%m/%Y")

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

    for filetype in export_mode:
        if filetype == "html":
            context = {
                "plots": plots,
                "export_mode": filetype,
                "datetime": completion_time,
            }
            jinja_render = template.render(context)
            out = output + ".html"
            with open(out, mode="w", encoding="utf-8") as f:
                f.write(jinja_render)
        else:
            context = {
                "plots": plots,
                "export_mode": filetype,
                "datetime": completion_time,
            }
            jinja_render = template.render(context)
            out = output + ".pdf"
            convert_html_to_pdf(jinja_render, out)

    print(f"Done! Your report can be found in {out}")

"""
Code in this script is used to generate the HTML output report
The functions are called by the main script RibosomeProfiler.py if the user specifies the --html flag
"""

from file_parser import parse_bam
from modules import *
from plots import *
import plotly.io as pio
from jinja2 import Environment, PackageLoader

env = Environment(loader=PackageLoader("RibosomeProfiler"), autoescape=False)
out = "RibosomeProfiler_report.html"
bam_file = input("Please provide the path to the bam file.\n")

read_df = parse_bam(bam_file, 10000)
read_length_dict = read_length_distribution(read_df)
read_length_fig = plot_read_length_distribution(read_length_dict, dict())
ligation_bias_dict = ligation_bias_distribution(read_df)
ligation_bias_fig = plot_ligation_bias_distribution(ligation_bias_dict, dict())
nucleotide_composition_dict = nucleotide_composition(read_df)
nucleotide_composition_fig = plot_nucleotide_composition(
    nucleotide_composition_dict, dict()
)

plots = [
    {
        "name": "Read Length Distribution",
        "description": "A plot showcasing the read length distribution of the reads",
        "fig": pio.to_html(read_length_fig, full_html=False),
    },
    {
        "name": "Ligation Bias Distribution",
        "description": "A plot showcasing the ligation bias distribution of the reads",
        "fig": pio.to_html(ligation_bias_fig, full_html=False),
    },
    {
        "name": "Nucleotide Composition",
        "description": "A plot showcasing the nucleotide composition of the reads",
        "fig": pio.to_html(nucleotide_composition_fig, full_html=False),
    },
]

template = env.get_template("base.html")
context = {"plots": plots}
with open(out, mode="w", encoding="utf-8") as f:
    f.write(template.render(context))

print(f"Done! Your report can be found in {out}")

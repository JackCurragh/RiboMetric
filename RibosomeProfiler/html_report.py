"""
Code in this script is used to generate the HTML output report
The functions are called by the main script RibosomeProfiler.py if the user specifies the --html flag

NOTE: This is a draft with the goal of creating the html report through this script.
"""

from file_parser import parse_bam
from modules import *
from plots import *
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

env = Environment(loader=FileSystemLoader("RibosomeProfiler/templates"),autoescape=False)
out = "RibosomeProfiler_report.html"
bam_file = input("Please provide the path to the bam file.\n")

read_df_pre = parse_bam(bam_file, 10000)
read_df = read_df_pre.loc[read_df_pre.index.repeat(read_df_pre['count'])].reset_index(drop=True)
read_length_dict = read_length_distribution(read_df)
plot_read_length_dict = plot_read_length_distribution(read_length_dict, dict())
ligation_bias_dict = ligation_bias_distribution(read_df)
plot_ligation_bias_dict = plot_ligation_bias_distribution(ligation_bias_dict, dict())
nucleotide_composition_dict = nucleotide_composition(read_df)
plot_nucleotide_composition_dict = plot_nucleotide_composition(
    nucleotide_composition_dict, dict()
)
plots = [plot_read_length_dict, plot_ligation_bias_dict, plot_nucleotide_composition_dict]

completion_time = datetime.now().strftime("%H:%M:%S %d/%m/%Y")
#completion_time = "test"

template = env.get_template("base.html")
context = {"plots": plots, "datetime": completion_time}
with open(out, mode="w", encoding="utf-8") as f:
    f.write(template.render(context))

print(f"Done! Your report can be found in {out}")

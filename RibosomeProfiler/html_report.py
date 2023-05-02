"""
Code in this script is used to generate the HTML output report

The functions are called by the main script RibosomeProfiler.py if the user
specifies the --html flag
"""

import pandas as pd
import plotly
from modules import read_length_distribution, ligation_bias_distribution


from plots import plot_read_length_distribution
from plots import plot_ligation_bias_distribution

out = "RibosomeProfiler-report.html"  # temp


def create_html_report(
    read_df: pd.DataFrame,
):  # add variables such as outdir, filename?
    """
    create a html report of RibosomeProfiler containing QC plots and more

    Inputs:
        read_df: Dataframe containing the read information
        out: Path where the html file should be stored (Default = .)

    Outputs:
        html: A html page containing the QC report from RibosomeProfiler
    """
    # create HTML page
    html_page = str()
    head = "<header>"
    head += "<h1>Ribosome Profiler Report</h1>"
    head += "</header>"
    html_page += head

    # Read Length Distribution
    body = "<body>"
    body += "<h2>Read Length Distribution</h2>"
    body += "<p>This shows how the read lengths are distributed.</p>"

    # create plot
    read_length_dict = read_length_distribution(read_df)
    read_length_plot = plot_read_length_distribution(read_length_dict, dict())
    body += plotly.io.to_html(read_length_plot, full_html=False)

    # Ligation Bias
    body += "<h2>Ligation Bias</h2>"
    body += """<p>This shows the distribution of the first two nucleotides
    in the reads to check for ligation bias.</p>"""

    # create plot
    ligation_bias_dict = ligation_bias_distribution(read_df)
    ligation_bias_plot = plot_ligation_bias_distribution(ligation_bias_dict,
                                                         dict()
                                                         )
    body += plotly.io.to_html(ligation_bias_plot, full_html=False)

    body += "</body>"
    html_page += body

    # save HTML page to file
    with open(out, "w") as f:  # temp outdir
        f.write(str(html_page))

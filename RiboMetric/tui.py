"""
RiboMetric TUI Viewer

A comprehensive Text User Interface for exploring RiboMetric analysis results.
Supports viewing JSON output files with interactive navigation and visualization.
"""

import json
import sys
from pathlib import Path
from typing import Dict, Any, Optional
import base64
import tempfile

from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Static, TabbedContent, TabPane
try:
    from textual.widgets import Image as TUIImage  # textual >=0.25
except Exception:  # pragma: no cover
    TUIImage = None
from textual.binding import Binding
from textual.screen import Screen

from .plots import (
    plot_read_frame_distribution,
    plot_metagene_heatmap,
    plot_read_length_distribution,
    plot_terminal_nucleotide_bias_distribution,
)


class RiboMetricData:
    """Wrapper class for RiboMetric JSON data with convenient accessors"""

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.data = self._load_data()
        self.results = self.data.get("results", {})
        self.config = self.data.get("config", {})
        self.metrics = self.results.get("metrics", {})
        self.mode = self.results.get("mode", "unknown")

    def _load_data(self) -> dict:
        """Load JSON data from file"""
        with open(self.filepath, 'r') as f:
            return json.load(f)

    def get_sample_name(self) -> str:
        """Extract sample name from filename"""
        return self.filepath.stem.replace("_RiboMetric", "").replace("_data", "")

    def get_total_reads(self) -> int:
        """Calculate total number of reads"""
        return sum(self.results.get("read_length_distribution", {}).values())

    def get_global_metric(self, metric_name: str) -> Optional[float]:
        """Get global value for a metric"""
        metric = self.metrics.get(metric_name)
        if isinstance(metric, dict):
            return metric.get("global")
        return metric

    def get_qc_status(self) -> str:
        """Determine overall QC status based on key metrics"""
        thresholds = {
            "periodicity_dominance": 0.5,
            "prop_reads_CDS": 0.5,
        }

        failures = 0
        for metric_name, threshold in thresholds.items():
            value = self.get_global_metric(metric_name)
            if value is not None and value < threshold:
                failures += 1

        if failures == 0:
            return "PASS"
        elif failures <= 1:
            return "WARNING"
        else:
            return "FAIL"


def render_summary(data: RiboMetricData) -> str:
    """Render summary view content"""
    lines = []

    # QC Status
    qc_status = data.get_qc_status()
    status_color = {"PASS": "green", "WARNING": "yellow", "FAIL": "red"}[qc_status]
    lines.append(f"[bold]QC Status:[/bold] [{status_color}]{qc_status}[/{status_color}]\n")

    # Sample info
    lines.append("[bold cyan]Sample Information:[/bold cyan]")
    lines.append(f"  Sample: {data.get_sample_name()}")
    lines.append(f"  Mode: {data.mode}")
    lines.append(f"  Total Reads: {data.get_total_reads():,}\n")

    # Read length stats
    rl_dist = data.results.get("read_length_distribution", {})
    if rl_dist:
        lengths = sorted([int(k) for k in rl_dist.keys()])
        most_common = max(rl_dist.items(), key=lambda x: x[1])
        lines.append(f"  Read Length Range: {min(lengths)}-{max(lengths)} nt")
        lines.append(f"  Most Common Length: {most_common[0]} nt ({most_common[1]:,} reads)\n")

    # Key metrics
    lines.append("[bold cyan]Key Metrics:[/bold cyan]")

    key_metrics = [
        ("periodicity_dominance", "Periodicity Dominance"),
        ("prop_reads_CDS", "CDS Proportion"),
        ("CDS_coverage_metric", "CDS Coverage"),
        ("read_length_distribution_IQR_metric", "Read Length IQR"),
    ]

    for metric_key, display_name in key_metrics:
        value = data.get_global_metric(metric_key)
        if value is not None:
            bar_width = 30
            filled = int(value * bar_width)
            bar = "█" * filled + "░" * (bar_width - filled)
            lines.append(f"  {display_name:30s} {value:6.3f}  {bar}")

    return "\n".join(lines)


def render_read_lengths(data: RiboMetricData) -> str:
    """Render read length distribution view"""
    lines = []
    lines.append("[bold cyan]Read Length Distribution:[/bold cyan]\n")

    rl_dist = data.results.get("read_length_distribution", {})
    if not rl_dist:
        return "No read length distribution data available"

    sorted_lengths = sorted([(int(k), v) for k, v in rl_dist.items()])
    total_reads = sum(v for _, v in sorted_lengths)
    max_count = max(v for _, v in sorted_lengths)

    lines.append(f"{'Length':>8s} {'Count':>12s} {'Proportion':>12s}  {'Distribution':40s}")
    lines.append("─" * 80)

    for length, count in sorted_lengths:
        proportion = count / total_reads
        bar_width = 35
        filled = int((count / max_count) * bar_width)
        bar = "█" * filled
        lines.append(f"{length:>8d} {count:>12,d} {proportion:>11.1%}  {bar}")

    # Per-length periodicity if available
    if "periodicity_dominance" in data.metrics:
        lines.append("\n[bold cyan]Per-Length Periodicity:[/bold cyan]\n")
        lines.append(f"{'Length':>8s} {'Reads':>12s} {'Periodicity':>12s}  {'Score':20s}")
        lines.append("─" * 60)

        periodicity = data.metrics["periodicity_dominance"]
        for length, count in sorted_lengths[:15]:  # Show first 15
            length_str = str(length)
            period_val = periodicity.get(length_str, 0)
            if isinstance(period_val, (int, float)):
                bar_width = 15
                filled = int(period_val * bar_width)
                bar = "█" * filled + "░" * (bar_width - filled)
                lines.append(f"{length:>8d} {count:>12,d} {period_val:>11.3f}  {bar}")

    return "\n".join(lines)


def render_periodicity(data: RiboMetricData) -> str:
    """Render periodicity analysis view"""
    lines = []

    if "periodicity_dominance" not in data.metrics:
        return "Periodicity data not available (annotation mode required)"

    periodicity = data.metrics["periodicity_dominance"]
    global_score = periodicity.get("global", 0)

    score_color = "green" if global_score > 0.7 else "yellow" if global_score > 0.5 else "red"
    lines.append(f"[bold]Global Periodicity Score:[/bold] [{score_color}]{global_score:.3f}[/{score_color}]\n")

    # Per-length breakdown
    per_length = {k: v for k, v in periodicity.items() if k != "global" and isinstance(v, (int, float))}
    if per_length:
        lines.append("[bold cyan]Per-Length Periodicity:[/bold cyan]\n")
        lines.append(f"{'Length':>8s} {'Score':>10s}  {'Visualization':30s}  {'Status':10s}")
        lines.append("─" * 70)

        sorted_lengths = sorted([(int(k), v) for k, v in per_length.items()])
        for length, score in sorted_lengths:
            bar_width = 25
            filled = int(score * bar_width)
            bar = "█" * filled + "░" * (bar_width - filled)

            if score > 0.7:
                status = "[green]✓ Good[/green]"
            elif score > 0.5:
                status = "[yellow]⚠ Weak[/yellow]"
            else:
                status = "[red]✗ Poor[/red]"

            lines.append(f"{length:>8d} {score:>10.3f}  {bar}  {status}")

    return "\n".join(lines)


def render_regions(data: RiboMetricData) -> str:
    """Render regional distribution view"""
    lines = []
    lines.append("[bold cyan]Regional Read Distribution:[/bold cyan]\n")

    prop_cds = data.get_global_metric("prop_reads_CDS")
    prop_leader = data.get_global_metric("prop_reads_leader")
    prop_trailer = data.get_global_metric("prop_reads_trailer")

    if prop_cds is None:
        return "Regional distribution data not available"

    lines.append(f"{'Region':>15s} {'Proportion':>12s}  {'Distribution':40s}")
    lines.append("─" * 75)

    regions = [
        ("5' Leader", prop_leader),
        ("CDS", prop_cds),
        ("3' Trailer", prop_trailer)
    ]

    for region_name, proportion in regions:
        if proportion is not None:
            bar_width = 35
            filled = int(proportion * bar_width)
            bar = "█" * filled + "░" * (bar_width - filled)
            lines.append(f"{region_name:>15s} {proportion:>11.1%}  {bar}")

    # Ratios
    ratio_cds_leader = data.get_global_metric("ratio_cds:leader")
    ratio_cds_trailer = data.get_global_metric("ratio_cds:trailer")

    if ratio_cds_leader is not None or ratio_cds_trailer is not None:
        lines.append("\n[bold cyan]Regional Ratios:[/bold cyan]\n")
        if ratio_cds_leader is not None:
            lines.append(f"  CDS : 5' Leader  = {ratio_cds_leader:.2f}")
        if ratio_cds_trailer is not None:
            lines.append(f"  CDS : 3' Trailer = {ratio_cds_trailer:.2f}")

    return "\n".join(lines)


def render_all_metrics(data: RiboMetricData) -> str:
    """Render all metrics table"""
    lines = []
    lines.append("[bold cyan]All Calculated Metrics:[/bold cyan]\n")
    lines.append(f"{'Metric':50s} {'Type':12s} {'Value':>15s}")
    lines.append("─" * 85)

    for metric_name, metric_value in sorted(data.metrics.items()):
        if isinstance(metric_value, dict):
            # Show global value
            if "global" in metric_value:
                lines.append(f"{metric_name:50s} {'Global':12s} {metric_value['global']:>15.4f}")
            # Show a few per-length values
            per_length = {k: v for k, v in metric_value.items() if k != "global"}
            count = 0
            for k, v in sorted(per_length.items(), key=lambda x: str(x[0]))[:3]:
                display_name = f"  {metric_name} [RL {k}]" if count == 0 else f"    [RL {k}]"
                lines.append(f"{display_name:50s} {'Per-length':12s} {v:>15.4f}")
                count += 1
            if len(per_length) > 3:
                lines.append(f"{'    ... and ' + str(len(per_length) - 3) + ' more':50s} {''}")
        elif isinstance(metric_value, (int, float)):
            val_str = f"{metric_value:.4f}" if isinstance(metric_value, float) else str(metric_value)
            lines.append(f"{metric_name:50s} {'Scalar':12s} {val_str:>15s}")

    return "\n".join(lines)


def render_raw_data(data: RiboMetricData) -> str:
    """Render raw JSON data"""
    json_str = json.dumps(data.data, indent=2)
    # Limit to first 100 lines for display
    lines = json_str.split('\n')[:100]
    if len(json_str.split('\n')) > 100:
        lines.append("... (truncated, use text editor to view full file)")
    return "\n".join(lines)


class HelpScreen(Screen):
    """Help screen showing keyboard shortcuts"""

    BINDINGS = [
        ("escape", "app.pop_screen", "Close"),
        ("q", "app.pop_screen", "Close"),
    ]

    def compose(self) -> ComposeResult:
        help_text = """[bold cyan]RiboMetric TUI Viewer - Help[/bold cyan]

[bold]Navigation:[/bold]
  [cyan]↑/↓[/cyan]        Navigate through content
  [cyan]←/→[/cyan]        Switch between tabs
  [cyan]Tab[/cyan]         Cycle through interactive elements
  [cyan]1-6[/cyan]         Quick jump to specific views

[bold]Views:[/bold]
  [cyan]1[/cyan]  Summary      - Overall QC status and key metrics
  [cyan]2[/cyan]  Read Lengths - Read length distribution and statistics
  [cyan]3[/cyan]  Periodicity  - Frame preference and periodicity analysis
  [cyan]4[/cyan]  Regions      - Regional read distribution (CDS, UTRs)
  [cyan]5[/cyan]  Metrics      - Complete table of all calculated metrics
  [cyan]6[/cyan]  Raw Data     - Raw JSON data view

[bold]Commands:[/bold]
  [cyan]h, ?[/cyan]       Show this help screen
  [cyan]q[/cyan]          Quit the application
  [cyan]r[/cyan]          Reload data from file

Press [cyan]Esc[/cyan] or [cyan]q[/cyan] to close this help screen.
"""
        yield Static(help_text)


class RiboMetricTUI(App):
    """Main TUI application for RiboMetric viewer"""

    CSS = """
    Static {
        padding: 1;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit", priority=True),
        Binding("h", "help", "Help"),
        Binding("question_mark", "help", "Help"),
        Binding("r", "reload", "Reload"),
        Binding("1", "switch_tab('summary')", "Summary"),
        Binding("2", "switch_tab('reads')", "Reads"),
        Binding("3", "switch_tab('periodicity')", "Periodicity"),
        Binding("4", "switch_tab('regions')", "Regions"),
        Binding("5", "switch_tab('metrics')", "Metrics"),
        Binding("6", "switch_tab('raw')", "Raw Data"),
    ]

    TITLE = "RiboMetric TUI Viewer"

    def __init__(self, json_file: str, **kwargs):
        super().__init__(**kwargs)
        self.json_file = json_file
        # Load data immediately
        try:
            self.data = RiboMetricData(self.json_file)
        except Exception as e:
            print(f"Error loading file: {e}")
            sys.exit(1)
        # Temp dir and plot cache
        self._plot_tmpdir = Path(tempfile.mkdtemp(prefix="ribometric_tui_"))
        self._plots: Dict[str, Optional[Path]] = {}

    def _write_image_from_b64(self, b64: str, name: str) -> Optional[Path]:
        if not b64:
            return None
        try:
            img_bytes = base64.b64decode(b64)
            path = self._plot_tmpdir / f"{name}.jpg"
            path.write_bytes(img_bytes)
            return path
        except Exception:
            return None

    def _build_plots(self) -> None:
        """Build plot images for available data; safe if disabled in CI."""
        try:
            cfg = self.data.config or {}
            res = self.data.results or {}
            self._plots.clear()
            if res.get("read_frame_distribution"):
                p = plot_read_frame_distribution(res["read_frame_distribution"], cfg)
                self._plots["frames_plot"] = self._write_image_from_b64(p.get("fig_image", ""), "frames_plot")
            if res.get("read_length_distribution"):
                p = plot_read_length_distribution(res["read_length_distribution"], cfg)
                self._plots["lengths_plot"] = self._write_image_from_b64(p.get("fig_image", ""), "lengths_plot")
            if res.get("terminal_nucleotide_bias_distribution"):
                p = plot_terminal_nucleotide_bias_distribution(res["terminal_nucleotide_bias_distribution"], cfg)
                self._plots["ligation_plot"] = self._write_image_from_b64(p.get("fig_image", ""), "ligation_plot")
            if res.get("metagene_profile"):
                p = plot_metagene_heatmap(res["metagene_profile"], cfg)
                self._plots["metagene_plot"] = self._write_image_from_b64(p.get("fig_image", ""), "metagene_plot")
        except Exception:
            # Non-fatal if images disabled or deps missing
            self._plots.clear()

    def on_mount(self) -> None:
        """Set subtitle when app starts"""
        self.sub_title = f"File: {Path(self.json_file).name}"
        # Build plots on start
        self._build_plots()

    def compose(self) -> ComposeResult:
        """Create the main layout"""
        yield Header()

        with TabbedContent(initial="summary"):
            with TabPane("Summary", id="summary"):
                yield Static(render_summary(self.data))

            with TabPane("Read Lengths", id="reads"):
                yield Static(render_read_lengths(self.data))

            with TabPane("Periodicity", id="periodicity"):
                yield Static(render_periodicity(self.data))

            with TabPane("Regions", id="regions"):
                yield Static(render_regions(self.data))

            with TabPane("All Metrics", id="metrics"):
                yield Static(render_all_metrics(self.data))

            with TabPane("Raw Data", id="raw"):
                yield Static(render_raw_data(self.data))

            # Optional visual plots (images)
            with TabPane("Lengths (Plot)", id="lengths_plot"):
                if TUIImage and self._plots.get("lengths_plot"):
                    yield TUIImage(self._plots["lengths_plot"])  # type: ignore[arg-type]
                else:
                    yield Static("Plot image not available.")

            with TabPane("Frames (Plot)", id="frames_plot"):
                if TUIImage and self._plots.get("frames_plot"):
                    yield TUIImage(self._plots["frames_plot"])  # type: ignore[arg-type]
                else:
                    yield Static("Plot image not available.")

            with TabPane("Lig. Bias (Plot)", id="ligation_plot"):
                if TUIImage and self._plots.get("ligation_plot"):
                    yield TUIImage(self._plots["ligation_plot"])  # type: ignore[arg-type]
                else:
                    yield Static("Plot image not available (requires sequence background).")

            with TabPane("Metagene (Plot)", id="metagene_plot"):
                if TUIImage and self._plots.get("metagene_plot"):
                    yield TUIImage(self._plots["metagene_plot"])  # type: ignore[arg-type]
                else:
                    yield Static("Plot image not available (requires annotation).")

        yield Footer()

    def action_help(self) -> None:
        """Show help screen"""
        self.push_screen(HelpScreen())

    def action_reload(self) -> None:
        """Reload data from file"""
        try:
            self.data = RiboMetricData(self.json_file)
            self.notify("Data reloaded successfully", severity="information")
            self._build_plots()

            # Update all tab contents
            tabs = self.query_one(TabbedContent)
            summary_static = tabs.get_pane("summary").query_one(Static)
            summary_static.update(render_summary(self.data))

            reads_static = tabs.get_pane("reads").query_one(Static)
            reads_static.update(render_read_lengths(self.data))

            periodicity_static = tabs.get_pane("periodicity").query_one(Static)
            periodicity_static.update(render_periodicity(self.data))

            regions_static = tabs.get_pane("regions").query_one(Static)
            regions_static.update(render_regions(self.data))

            metrics_static = tabs.get_pane("metrics").query_one(Static)
            metrics_static.update(render_all_metrics(self.data))

            raw_static = tabs.get_pane("raw").query_one(Static)
            raw_static.update(render_raw_data(self.data))

            # Update image panes
            if TUIImage:
                try:
                    lpane = tabs.get_pane("lengths_plot")
                    if self._plots.get("lengths_plot"):
                        lpane.remove_children()
                        lpane.mount(TUIImage(self._plots["lengths_plot"]))  # type: ignore[arg-type]
                except Exception:
                    pass
                try:
                    fpane = tabs.get_pane("frames_plot")
                    if self._plots.get("frames_plot"):
                        # Replace children with a fresh Image widget
                        fpane.remove_children()
                        fpane.mount(TUIImage(self._plots["frames_plot"]))  # type: ignore[arg-type]
                except Exception:
                    pass
                try:
                    ligpane = tabs.get_pane("ligation_plot")
                    if self._plots.get("ligation_plot"):
                        ligpane.remove_children()
                        ligpane.mount(TUIImage(self._plots["ligation_plot"]))  # type: ignore[arg-type]
                except Exception:
                    pass
                try:
                    mpane = tabs.get_pane("metagene_plot")
                    if self._plots.get("metagene_plot"):
                        mpane.remove_children()
                        mpane.mount(TUIImage(self._plots["metagene_plot"]))  # type: ignore[arg-type]
                except Exception:
                    pass

        except Exception as e:
            self.notify(f"Error reloading data: {e}", severity="error")

    def action_switch_tab(self, tab_id: str) -> None:
        """Switch to specific tab"""
        tabs = self.query_one(TabbedContent)
        tabs.active = tab_id


def run_tui(json_file: str):
    """Launch the TUI application"""
    app = RiboMetricTUI(json_file)
    app.run()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python tui.py <RiboMetric_JSON_file>")
        sys.exit(1)

    run_tui(sys.argv[1])

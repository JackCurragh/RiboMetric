from RiboMetric.html_report import build_report_context


def test_build_report_context_groups_and_cards():
    summary = {
        "metrics": [
            {
                "key": "periodicity_dominance_score",
                "metric": "periodicity_dominance",
                "score": 0.82,
            },
            {"key": "cds_enrichment_score", "metric": "cds_enrichment_ratio", "score": 0.76},
            {"key": "fragment_uniqueness_score", "metric": "duplicate_rate", "score": 0.12},
            {"key": "rpf_unique_mapping_score", "metric": "rpf_multimapper_rate", "score": 0.34},
            {
                "key": "terminal_evenness_kl_5prime_score",
                "metric": "terminal_bias_kl_5prime",
                "score": 0.61,
            },
        ]
    }

    context = build_report_context(summary)

    # 0.12 is a failing score full stop. Under the old naming this row was
    # "duplicate_rate" and had to be special-cased as lower-is-better before
    # its status could be read; a score key needs no such caveat.
    assert context["overall_status"] == "fail"
    # Cards are addressed by score key; grouping follows the raw metric.
    assert [card["key"] for card in context["top_cards"]] == [
        "periodicity_dominance_score",
        "cds_enrichment_score",
        "fragment_uniqueness_score",
        "rpf_unique_mapping_score",
    ]
    groups = {group["name"]: group for group in context["grouped_metrics"]}
    assert "Mapping" in groups
    assert "Periodicity" in groups
    assert "Codon / Sequence" in groups


def test_every_score_is_labelled_higher_is_better():
    """Scores are the only 0-1 numbers in the report and they all point the
    same way, so the report never has to caveat a direction."""
    summary = {
        "metrics": [
            {
                "key": "fragment_uniqueness_score",
                "metric": "duplicate_rate",
                "score": 0.8,
                "raw": 0.2,
                "status": "PASS",
                "tier": 3,
            },
            {
                "key": "periodicity_dominance_score",
                "metric": "periodicity_dominance",
                "score": 0.8,
                "raw": 0.8,
                "status": "PASS",
                "tier": 1,
            },
        ]
    }

    context = build_report_context(summary)
    rows = [metric for group in context["grouped_metrics"] for metric in group["metrics"]]
    assert rows, "expected the legacy grouped view to be populated"
    for row in rows:
        assert row["direction"] == "Higher is better"
        assert row["status"] == "pass"

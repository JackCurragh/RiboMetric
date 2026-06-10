from RiboMetric.html_report import build_report_context


def test_build_report_context_groups_and_cards():
    summary = {
        "metrics": [
            {"key": "periodicity_dominance", "name": "Periodicity dominance", "score": 0.82},
            {"key": "prop_reads_CDS", "name": "Prop reads CDS", "score": 0.76},
            {"key": "duplicate_rate", "name": "Duplicate rate", "score": 0.12},
            {"key": "multimapper_rate", "name": "Multimapper rate", "score": 0.34},
            {"key": "terminal_bias_kl_5prime_score", "name": "Bias", "score": 0.61},
        ]
    }

    context = build_report_context(summary)

    assert context["overall_status"] == "warn"
    assert [card["key"] for card in context["top_cards"]] == [
        "periodicity_dominance",
        "prop_reads_CDS",
        "duplicate_rate",
        "multimapper_rate",
    ]
    groups = {group["name"]: group for group in context["grouped_metrics"]}
    assert "Mapping" in groups
    assert "Periodicity" in groups
    assert "Codon / Sequence" in groups


def test_build_report_context_respects_lower_is_better_status():
    summary = {
        "metrics": [
            {"key": "duplicate_rate", "name": "Duplicate rate", "score": 0.8},
            {"key": "periodicity_dominance", "name": "Periodicity dominance", "score": 0.8},
        ]
    }

    context = build_report_context(summary)
    statuses = {
        metric["key"]: metric["status"]
        for group in context["grouped_metrics"]
        for metric in group["metrics"]
    }

    assert statuses["duplicate_rate"] == "fail"
    assert statuses["periodicity_dominance"] == "pass"

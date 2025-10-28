# Testing Guide for RiboMetric

This guide explains how to run and write tests for RiboMetric.

## Quick Start

### Install Test Dependencies

```bash
pip install -r requirements_test.txt
```

### Run All Tests

```bash
pytest
```

### Run Tests with Coverage

```bash
pytest --cov=RiboMetric --cov-report=term-missing
```

### Run Specific Test Files

```bash
# Test the new metric selection feature
pytest tests/test_metric_selection.py -v

# Test metrics
pytest tests/test_metrics_improved.py -v

# Test modules
pytest tests/test_modules.py -v
```

### Run Tests by Marker

```bash
# Run only unit tests (fast)
pytest -m unit

# Run only integration tests
pytest -m integration

# Skip slow tests
pytest -m "not slow"
```

## Test Structure

```
tests/
├── conftest.py                    # Pytest fixtures and configuration
├── test_RiboMetric.py            # Integration tests for main workflows
├── test_metric_selection.py       # Tests for default/optional metric selection
├── test_metrics.py               # Original metric tests
├── test_metrics_improved.py       # Enhanced metric tests with fixtures
├── test_modules.py               # Tests for module functions
├── test_plots.py                 # Tests for plotting functions
├── test_parser.py                # Tests for file parsing
└── test_data/                    # Test data files
    ├── test.bam
    ├── test.csv
    ├── test_annotation.tsv
    └── ...
```

## Writing Tests

### Use Fixtures

Fixtures are defined in `conftest.py` and can be used in any test:

```python
def test_my_metric(sample_read_length_dict):
    """Use the sample_read_length_dict fixture"""
    result = my_metric_function(sample_read_length_dict)
    assert result > 0
```

Available fixtures:
- `test_data_dir` - Path to test data directory
- `config_path` - Path to config file
- `sample_read_df` - Sample read dataframe
- `sample_annotation_df` - Sample annotation dataframe
- `sample_read_length_dict` - Sample read length distribution
- `sample_read_frame_dict` - Sample read frame distribution
- `sample_metagene_profile` - Sample metagene profile
- `sample_sequence_background` - Sample sequence background frequencies
- `default_config` - Default configuration dictionary

### Test Organization

Organize tests into classes by functionality:

```python
class TestMetricName:
    """Tests for specific metric"""

    def test_normal_case(self):
        """Test normal operation"""
        pass

    def test_edge_case(self):
        """Test edge cases"""
        pass

    def test_error_handling(self):
        """Test error conditions"""
        pass
```

### Parameterized Tests

Use `@pytest.mark.parametrize` for testing multiple inputs:

```python
@pytest.mark.parametrize("input_val,expected", [
    (10, 100),
    (20, 400),
    (30, 900),
])
def test_square(input_val, expected):
    assert input_val ** 2 == expected
```

### Mark Tests

Use markers to categorize tests:

```python
@pytest.mark.slow
def test_full_pipeline():
    """Long-running integration test"""
    pass

@pytest.mark.unit
def test_single_function():
    """Fast unit test"""
    pass

@pytest.mark.requires_data
def test_with_real_data(test_data_dir):
    """Test that needs real data files"""
    pass
```

## Test Coverage

### Generate Coverage Report

```bash
pytest --cov=RiboMetric --cov-report=html
```

Then open `htmlcov/index.html` in your browser.

### Coverage Goals

- **Overall**: Aim for >80% coverage
- **Critical modules**: >90% coverage
  - `metrics.py`
  - `modules.py`
  - `qc.py`
- **UI/CLI**: >70% coverage
  - `cli.py`
  - `arg_parser.py`

## Continuous Integration

Tests run automatically on GitHub Actions for:
- Python 3.8, 3.9, 3.10
- On pushes to main, dev branches
- On pull requests

View results at: https://github.com/JackCurragh/RiboMetric/actions

## Common Test Scenarios

### Testing Metric Calculations

```python
def test_metric_range():
    """Ensure metric is in valid range"""
    result = my_metric(sample_data)
    assert 0 <= result <= 1

def test_metric_reproducibility():
    """Ensure metric gives same result"""
    result1 = my_metric(sample_data)
    result2 = my_metric(sample_data)
    assert result1 == result2

def test_metric_edge_cases():
    """Test with empty/zero/extreme inputs"""
    assert my_metric({}) == 0
    assert my_metric({0: 0}) == 0
```

### Testing Default vs Optional Metrics

```python
def test_default_metrics_only(default_config):
    """Test that optional metrics are skipped"""
    config = default_config.copy()
    config["enabled_metrics"] = config["metrics"]["default"]

    assert should_calculate_metric("periodicity_dominance", config)
    assert not should_calculate_metric("periodicity_fourier", config)
```

### Testing File I/O

```python
def test_file_parsing(test_data_dir):
    """Test parsing of data files"""
    file_path = test_data_dir / "test.csv"
    result = parse_file(file_path)
    assert not result.empty
    assert "column_name" in result.columns
```

## Debugging Tests

### Run with verbose output

```bash
pytest -vv
```

### Run with print statements

```bash
pytest -s
```

### Run single test

```bash
pytest tests/test_file.py::TestClass::test_method -v
```

### Drop into debugger on failure

```bash
pytest --pdb
```

## Best Practices

1. **One assertion per test** (when possible) - Makes failures easier to diagnose
2. **Use descriptive test names** - `test_metric_handles_zero_division` not `test_1`
3. **Test edge cases** - Empty inputs, zero values, very large values
4. **Use fixtures** - Don't duplicate test data setup
5. **Keep tests fast** - Mock slow operations, mark slow tests
6. **Test behavior, not implementation** - Test what the function does, not how

## Troubleshooting

### ImportError: No module named 'RiboMetric'

Install the package in development mode:

```bash
pip install -e .
```

### Tests can't find test data

Use the `test_data_dir` fixture instead of hardcoded paths.

### Samtools not found

Install samtools:

```bash
# macOS
brew install samtools

# Ubuntu/Debian
sudo apt-get install samtools

# conda
conda install -c bioconda samtools
```

## Contributing

When adding new features:

1. Write tests first (TDD approach)
2. Ensure all tests pass
3. Maintain or improve coverage
4. Add docstrings to test functions
5. Update this guide if needed

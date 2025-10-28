#!/bin/bash
# Convenience script for running RiboMetric tests

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}RiboMetric Test Runner${NC}"
echo "======================"
echo ""

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest not found${NC}"
    echo "Install test dependencies with:"
    echo "  pip install -r requirements_test.txt"
    exit 1
fi

# Parse command line arguments
case "${1:-all}" in
    all)
        echo -e "${GREEN}Running all tests...${NC}"
        pytest -v
        ;;

    coverage)
        echo -e "${GREEN}Running tests with coverage...${NC}"
        pytest --cov=RiboMetric --cov-report=term-missing --cov-report=html
        echo ""
        echo -e "${GREEN}Coverage report saved to htmlcov/index.html${NC}"
        ;;

    fast)
        echo -e "${GREEN}Running fast tests only...${NC}"
        pytest -m "not slow" -v
        ;;

    unit)
        echo -e "${GREEN}Running unit tests...${NC}"
        pytest -m unit -v
        ;;

    integration)
        echo -e "${GREEN}Running integration tests...${NC}"
        pytest -m integration -v
        ;;

    metrics)
        echo -e "${GREEN}Running metric tests...${NC}"
        pytest tests/test_metrics.py tests/test_metrics_improved.py -v
        ;;

    selection)
        echo -e "${GREEN}Running metric selection tests...${NC}"
        pytest tests/test_metric_selection.py -v
        ;;

    modules)
        echo -e "${GREEN}Running module tests...${NC}"
        pytest tests/test_modules.py -v
        ;;

    parallel)
        echo -e "${GREEN}Running tests in parallel...${NC}"
        pytest -n auto -v
        ;;

    watch)
        echo -e "${GREEN}Running tests in watch mode...${NC}"
        pytest-watch
        ;;

    help)
        echo "Usage: ./run_tests.sh [command]"
        echo ""
        echo "Commands:"
        echo "  all         - Run all tests (default)"
        echo "  coverage    - Run tests with coverage report"
        echo "  fast        - Run fast tests only (skip slow tests)"
        echo "  unit        - Run unit tests only"
        echo "  integration - Run integration tests only"
        echo "  metrics     - Run metric tests"
        echo "  selection   - Run metric selection tests"
        echo "  modules     - Run module tests"
        echo "  parallel    - Run tests in parallel (requires pytest-xdist)"
        echo "  watch       - Run tests in watch mode (requires pytest-watch)"
        echo "  help        - Show this help message"
        echo ""
        echo "Examples:"
        echo "  ./run_tests.sh                  # Run all tests"
        echo "  ./run_tests.sh coverage         # Run with coverage"
        echo "  ./run_tests.sh fast             # Skip slow tests"
        ;;

    *)
        echo -e "${RED}Unknown command: $1${NC}"
        echo "Run './run_tests.sh help' for usage information"
        exit 1
        ;;
esac

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Tests completed successfully${NC}"
else
    echo ""
    echo -e "${RED}✗ Tests failed${NC}"
fi

exit $exit_code

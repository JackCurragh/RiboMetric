"""Console script for RibosomeProfiler."""
import sys

from RibosomeProfiler import argumnet_parser as p, main as m


def main():
    """Console script for RibosomeProfiler."""
    parser = p()
    args = parser.parse_args()

    m(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

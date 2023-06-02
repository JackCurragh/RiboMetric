"""Console script for RiboMetric."""
import sys

from .RiboMetric import argument_parser as p, main as m


def main():
    """Console script for RiboMetric."""
    parser = p()
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()

    m(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

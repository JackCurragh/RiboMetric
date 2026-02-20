"""Console script for RiboMetric."""
import sys
import logging
import argparse

from .RiboMetric import argument_parser as p, main as m
from . import __version__


def main() -> int:
    """Console script for RiboMetric."""
    parser: argparse.ArgumentParser = p()
    # Global flags
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase verbosity (-v, -vv)")
    parser.add_argument("-q", "--quiet", action="store_true", help="Reduce output")
    args = parser.parse_args()
    if getattr(args, "version", False):
        print(__version__)
        return 0

    # Configure logging
    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    if args.quiet:
        level = logging.ERROR
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")

    if not vars(args):
        parser.print_help()

    m(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

#!/usr/bin/env python

"""
Parse command-line arguments and external files.
"""

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments when executed with module options.

    Returns:
        argparse.Namespace: parse results.
    """
    arg_parser = argparse.ArgumentParser(
        description="A Python implementation of a part of compareDEtools."
    )
    arg_parser.add_argument(
        "-c",
        "--condition",
        nargs="*",
        type=Path,
        help="file path of simulation conditions",
    )
    args = arg_parser.parse_args()
    return args

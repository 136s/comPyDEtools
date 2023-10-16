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
        "condition",
        nargs="*",
        type=Path,
        help="file path of simulation conditions",
    )
    arg_parser.add_argument(
        "-g",
        "--generate",
        action="store_true",
        help="just generate simulated datasets",
    )
    arg_parser.add_argument(
        "-d",
        "--draw",
        action="store_true",
        help="When DE Analysis has done, just draw figures",
    )
    args = arg_parser.parse_args()
    return args

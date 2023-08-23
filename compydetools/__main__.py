#!/usr/bin/env python


from . import parser
from .condition import set_condition

if __name__ == "__main__":
    args = parser.parse_args()
    if args.condition:
        for cnd in args.condition:
            set_condition(cnd)

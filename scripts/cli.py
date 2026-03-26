import argparse


def parse_range(s: str):
    try:
        start, end = map(int, s.split('-'))
        return start, end
    except Exception:
        raise argparse.ArgumentTypeError("Range must be in START-END format (e.g. 5-10)")

#!/usr/bin/env python3
import sys

keep_from = int(sys.argv[2]) - 1
name = None
seq = []

with open(sys.argv[1]) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                s = "".join(seq)
                if "hxb2" not in name.lower():
                    s = "-" * keep_from + s[keep_from:]
                print(name)
                print(s)
            name = line
            seq = []
        else:
            seq.append(line.strip())

if name is not None:
    s = "".join(seq)
    if "hxb2" not in name.lower():
        s = "-" * keep_from + s[keep_from:]
    print(name)
    print(s)

#!/usr/bin/env python
"""Call library strandedness from RSeQC infer_experiment.py output."""
import argparse
import re
import sys

p = argparse.ArgumentParser()
p.add_argument("--infer", required=True)
p.add_argument("--mode", default="auto", help="auto|unstranded|forward|reverse")
p.add_argument("--threshold", type=float, default=0.8)
p.add_argument("--out", required=True)
a = p.parse_args()

valid = {"unstranded", "forward", "reverse"}
if a.mode != "auto":
    if a.mode not in valid:
        sys.exit(f"strandedness must be auto/unstranded/forward/reverse, got '{a.mode}'")
    call = a.mode
    print(f"Strandedness fixed by config: {call}")
else:
    text = open(a.infer).read()
    # First fraction = sense ("1++,1--,2+-,2-+" PE or "++,--" SE) -> forward
    # Second fraction = antisense -> reverse
    fr = re.findall(r'explained by "([^"]+)":\s*([0-9.eE+-]+)', text)
    if len(fr) != 2:
        sys.exit(f"Could not parse infer_experiment output: {a.infer}")
    fwd, rev = float(fr[0][1]), float(fr[1][1])
    tot = fwd + rev
    if tot == 0:
        sys.exit("infer_experiment reported no strand-informative reads")
    if fwd / tot >= a.threshold:
        call = "forward"
    elif rev / tot >= a.threshold:
        call = "reverse"
    else:
        call = "unstranded"
    print(f"forward={fwd:.4f} reverse={rev:.4f} -> {call}")

with open(a.out, "w") as fh:
    fh.write(call + "\n")

#!/usr/bin/env python3
"""make_variants.py - test variants of a real InfoSewer engine input file.

The real 125_21A model has only unpeakable loads and no patterns, so these
variants switch on the parts it never exercises. Each variant is run through the
real InfoSewer engine (ise_run) for the reference answer; swmm5-infosewer has to
match it. Changing a load's type or pattern is the only edit: geometry, load
values and options are left as they are.

  python3 -I make_variants.py <model.txt> <outdir>

  peak    every load peakable (type 1): Q = K * Qbase^P on accumulated flow
  cover   every load coverage (type 2) with units = 10 per load
  mixed   loads cycle through types 0, 1, 2 in file order
  diurnal every load unpeakable with a 24-value diurnal pattern DIURNAL (EPS)
  mixeddi mixed types, every load on DIURNAL (EPS and steady)
"""
import sys, os

DIURNAL = [0.55, 0.45, 0.40, 0.38, 0.40, 0.55, 0.85, 1.25, 1.45, 1.40, 1.30, 1.20,
           1.15, 1.10, 1.05, 1.05, 1.10, 1.20, 1.30, 1.35, 1.25, 1.05, 0.85, 0.70]


def rewrite(lines, mode):
    out, sect, k = [], None, 0
    for line in lines:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            sect = s
            out.append(line)
            if s == "[Patterns]" and mode in ("diurnal", "mixeddi"):
                for i in range(0, 24, 6):
                    out.append("DIURNAL 1\t" + "\t".join("%f" % v for v in DIURNAL[i:i + 6]) + "\n")
            continue
        if sect == "[Loading Manholes]" and s:
            t = line.rstrip("\n").split("\t")
            for g in range(10):
                i = 18 + 4 * g
                if i + 3 >= len(t) or t[i] == "-" or abs(float(t[i])) < 1e-12:
                    continue
                if mode == "peak":
                    t[i + 1] = "1"
                elif mode == "cover":
                    t[i + 1], t[i + 3] = "2", "10.000000"
                elif mode in ("mixed", "mixeddi"):
                    ty = k % 3
                    t[i + 1] = str(ty)
                    t[i + 3] = "10.000000" if ty == 2 else "0.000000"
                if mode in ("diurnal", "mixeddi"):
                    t[i + 2] = "DIURNAL"
                k += 1
            line = "\t".join(t) + "\n"
        out.append(line)
    return out


def main():
    src, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    lines = open(src, encoding="latin-1").readlines()
    base = os.path.splitext(os.path.basename(src))[0]
    for mode in ("peak", "cover", "mixed", "diurnal", "mixeddi"):
        p = os.path.join(outdir, "%s_%s.txt" % (base, mode))
        open(p, "w", encoding="latin-1").writelines(rewrite(lines, mode))
        print(p)


if __name__ == "__main__":
    main()

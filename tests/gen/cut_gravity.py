#!/usr/bin/env python3
"""cut_gravity.py - cut the gravity catchment upstream of one outlet out of an
InfoSewer engine input file, for gravity-only tests with real geometry, loads and
patterns. Wet wells, pumps, force mains, chambers, status and controls are dropped;
any pump or force main that discharged into the catchment is simply gone, so the
cut network is a test network, not the real system.

  python3 -I cut_gravity.py <model.txt> <outlet id> <out.txt>
"""
import sys, collections

DROP = ("Junction Chambers", "Wet Wells", "Force Mains", "Pumps", "Status", "Controls",
        "Flow Splits", "Extra Loads", "Load Node Head Type", "Initial Concentration")


def main():
    src, outlet, dst = sys.argv[1:4]
    lines = open(src, encoding="latin-1").read().split("\n")
    secs, order, cur = collections.OrderedDict(), [], None
    for l in lines:
        s = l.strip()
        if s.startswith("[") and s.endswith("]"):
            cur = s[1:-1]
            secs[cur] = []
            continue
        if cur is not None:
            secs[cur].append(l)
    up = collections.defaultdict(list)
    for l in secs.get("Pipes", []):
        t = l.split()
        if len(t) > 2:
            up[t[2]].append((t[0], t[1]))
    keep_nodes, keep_pipes, st = {outlet}, set(), [outlet]
    while st:
        n = st.pop()
        for pid, u in up[n]:
            keep_pipes.add(pid)
            if u not in keep_nodes:
                keep_nodes.add(u)
                st.append(u)
    out = []
    for name, body in secs.items():
        out.append("[%s]" % name)
        if name in DROP:
            out.append("")
            continue
        for l in body:
            t = l.split()
            if not t:
                out.append(l)
                continue
            if name == "Loading Manholes" and t[0] not in keep_nodes:
                continue
            if name == "Outlets" and t[0] != outlet:
                continue
            if name == "Outlet Heads" and t[0] != outlet:
                continue
            if name in ("Pipes", "Infiltrate Rates") and t[0] not in keep_pipes:
                continue
            out.append(l)
    open(dst, "w", encoding="latin-1").write("\n".join(out))
    print("%s: %d nodes, %d pipes -> %s" % (outlet, len(keep_nodes), len(keep_pipes), dst))


if __name__ == "__main__":
    main()

"""swmm_out.py - minimal EPA SWMM 5.x binary .out reader (link and node results).
Written by Claude (AI) for SWMM Enablement, 2026-11-09."""
import struct


def read(path):
    b = open(path, "rb").read()
    idpos, inpos, outpos, nper, err, magic = struct.unpack_from("<6i", b, len(b) - 24)
    _, ver, fu, ns, nn, nl, npol = struct.unpack_from("<7i", b, 0)
    p = idpos
    names = []
    for _ in range(ns + nn + nl + npol):
        n = struct.unpack_from("<i", b, p)[0]; p += 4
        names.append(b[p:p + n].decode("latin-1")); p += n
    nodes = names[ns:ns + nn]
    links = names[ns + nn:ns + nn + nl]
    p = inpos
    n = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * n + 4 * ns * n      # subcatch input props
    n = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * n + 4 * nn * n      # node input props
    n = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * n + 4 * nl * n      # link input props
    nsv = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * nsv
    nnv = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * nnv
    nlv = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * nlv
    nyv = struct.unpack_from("<i", b, p)[0]; p += 4 + 4 * nyv
    per = []
    rec = 8 + 4 * (ns * nsv + nn * nnv + nl * nlv + nyv)
    for k in range(nper):
        q = outpos + k * rec
        t = struct.unpack_from("<d", b, q)[0]; q += 8 + 4 * ns * nsv
        nv = struct.unpack_from("<%df" % (nn * nnv), b, q); q += 4 * nn * nnv
        lv = struct.unpack_from("<%df" % (nl * nlv), b, q)
        per.append(dict(time=t,
                        node={nodes[i]: nv[i * nnv:(i + 1) * nnv] for i in range(nn)},
                        link={links[i]: lv[i * nlv:(i + 1) * nlv] for i in range(nl)}))
    return dict(nodes=nodes, links=links, periods=per, flow_units=fu)

#!/usr/bin/env python2
"""
Generate sorting indices based one a good example granule.
It outputs C arrays.
"""

from __future__ import division, print_function, absolute_import

import numpy as np

SwathSize = 10

def writec(cols, widths):
    print("// DO NOT EDIT. This file was auto-generated by scripts/gensortindices.py")
    print("// using the latitude sorting indices of")
    print("// MOD03.A2015129.1540.005.2015131111937.hdf")

    print("")
    print("short SORT_WIDTHS[] = {")
    print("\t" + ", ".join(str(w) for w in widths))
    print("};")

    print("")
    print("short SORT_FIRST[][%d] = {" % (SwathSize,))
    for col in cols:
        print("\t{" + ", ".join(str(v) for v in col[:SwathSize]) + "},")
    print("};")

    print("")
    print("short SORT_MID[][%d] = {" % (SwathSize,))
    for col in cols:
        print("\t{" + ", ".join(str(v) for v in col[SwathSize:2*SwathSize]) + "},")
    print("};")

    print("")
    print("short SORT_LAST[][%d] = {" % (SwathSize,))
    for col in cols:
        print("\t{" + ", ".join(str(v) for v in col[-SwathSize:]) + "},")
    print("};")

def main():
    # sind.bin is a binary dump of the sorting indices
    sind = np.fromfile("sind.bin", dtype='int32').reshape((-1, 1354))
    coldiff = np.diff(sind, axis=1)
    spikes = np.mean(np.abs(coldiff), axis=0)

    nonzero = []
    for x in np.split(np.arange(len(spikes)), np.where(spikes == 0)[0]):
        if len(x) > 1:
            nonzero.append((x[1], x[-1], (x[1]+x[-1])//2))  # sic. ignore first one

    first = np.repeat(np.arange(sind.shape[0]//SwathSize)*SwathSize, SwathSize)
    prevend = 0
    prevmid = 0
    cols = []
    widths = []
    for start, end, mid in nonzero:
        col = np.median(sind[:,prevend:start], axis=1).astype('int32')
        cols.append(col - first)
        widths.append(mid-prevmid)
        prevend = end
        prevmid = mid

    col = np.median(sind[:,prevend:sind.shape[1]], axis=1).astype('int32')
    cols.append(col - first)
    widths.append(sind.shape[1]-prevmid)

    # check consistency of indices for each swath except first and last
    for col in cols:
        x = col[SwathSize:-SwathSize].reshape((-1, SwathSize))
        assert np.all(x[1:,:] == x[:-1,:])

    #newsind = []
    #for col, w in zip(cols, widths):
    #    newsind.append(np.tile(col[:,np.newaxis], w))
    #newsind = np.column_stack(newsind).astype('int32')
    #newsind.tofile("newsind.bin")

    writec(cols, widths)


if __name__ == '__main__':
    main()
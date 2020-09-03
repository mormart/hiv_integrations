import strutils, algorithm, tables

type
    Chromosomes* = enum
        chr1, chr2, chr3, chr4, chr5,
        chr6, chr7, chr8, chr9, chr10,
        chr11, chr12, chr13, chr14, chr15,
        chr16, chr17, chr18, chr19, chr20,
        chr21, chr22, chrX, chrY

type
    SqMatrix* = object
        width*: int
        data*: seq[int]

proc makeMatrix*(width: int): SqMatrix =
    result.width = width
    result.data = newSeq[int](width * width)

proc `inc`*(m: var SqMatrix, i, j: int) {.inline.} =
    assert i < m.width
    assert j < m.width
    m.data[i * m.width + j] += 1
    if i != j: m.data[j * m.width + i] += 1

proc getSeqlens*(filepath: string): OrderedTable[string, int] =
    let f = open(filepath)
    for line in f.lines:
        var fields = line.split("\t")
        result[fields[0]] = parseInt(fields[1])
    f.close()

proc findSites*(filepath: string, pattern: string, seqlens: OrderedTable): OrderedTable[Chromosomes, seq[int]] =
    var
        linenum = 0
        idx: int
        chrname: string
        ischr: bool
    let f = open(filepath)

    for line in f.lines:
        if line[0] == '>':
            if len(line) > 6 or line == ">chrM":
                ischr = false
                continue
            ischr = true
            if chrname != "": result[parseEnum[Chromosomes](chrname)].add(seqlens[chrname])
            chrname = line[1..line.high]
            result[parseEnum[Chromosomes](chrname)] = @[]
            linenum = 0
            continue       
        if ischr:
            idx = line.find(pattern)
            if idx != -1:
                result[parseEnum[Chromosomes](chrname)].add(linenum * 50 + idx)
            inc linenum

    if endOfFile(f): result[parseEnum[Chromosomes](chrname)].add(seqlens[chrname])
    f.close()

proc binSiteEnds*(fragments: OrderedTable, binsize: int): OrderedTable[Chromosomes, seq[int]] =
    var
        midpoint: int
        midpoints, ends: seq[int] = @[]

    for key in fragments.keys:
        ends = fragments[key]
        result[key] = @[]
        result[key].add(ends[^1])
        midpoints = @[]
        midpoints.add((ends[^1] + ends[^2]) div 2)

        for i in countdown(ends.high, 1):
            midpoint = (ends[i] + ends[i - 1]) div 2
            if midpoint div binsize == midpoints[^1] div binsize: continue
            midpoints.add(midpoint)
            result[key].add(ends[i])
        if (ends[0] div 2) div binsize != midpoints[^1] div binsize:
            result[key].add(ends[0])
        result[key] = result[key].reversed()
    result.sort(proc(x, y: tuple[key: Chromosomes, val: seq[int]]): int = cmp(ord(x.key), ord(y.key)))

proc binCumsum*(bins: OrderedTable): seq[int] =
    result.add(0)
    var sum = 0
    for key in bins.keys:
        inc(sum, len(bins[key]))
        result.add(sum)
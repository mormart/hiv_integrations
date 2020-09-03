import strutils, tables, makeFragments

proc findRegion(ends: seq[int], readpos, avgSize: int): int =
    var bestguess = readpos div avgSize

    if bestguess > ends.high:
        bestguess = ends.high
    if readpos < ends[bestguess]:
        while readpos < ends[bestguess]:
            if bestguess == 0: break
            dec bestguess
        inc bestguess
    else:
        while readpos > ends[bestguess]:
            inc bestguess
    return bestguess

echo "Getting seqlens..."
let seqlens = getSeqlens("hg38_seqlens.txt")

echo "Digesting genome..."
let fragments = findSites("hg38.fa", "GATC", seqlens = seqlens)

echo "Binning fragments..."
let bins = binSiteEnds(fragments, 100000)
let cumsum = binCumsum(bins)

echo "Finding contacts..."
let f = open("hic_jurkat2_allpairs.txt")

var
    readRegionIdx, mateRegionIdx: int
    contactMap = makeMatrix(cumsum[24])
    pairchr, matechr: Chromosomes
    counter = 0

for line in f.lines:
    if counter mod 10000000 == 0: echo counter
    inc counter
    var fields = line.split("\t")
    if fields[1].high > 4 or fields[3].high > 4 or fields[1] == "chrM" or fields[3] == "chrM": continue
    pairchr = parseEnum[Chromosomes](fields[1])
    if fields[3] == "=": matechr = pairchr
    else: matechr = parseEnum[Chromosomes](fields[3])
    # matechr = parseEnum[Chromosomes](fields[3])
    readRegionIdx = findRegion(bins[pairchr], parseInt(fields[2]), 100000) + cumsum[ord(pairchr)]
    mateRegionIdx = findRegion(bins[matechr], parseInt(fields[4]), 100000) + cumsum[ord(matechr)]
    inc(contactMap, readRegionIdx, mateRegionIdx)

f.close()

let output = open("jurkat_ALL_100K.txt", fmWrite)
for i in 0 .. contactMap.width - 1:
    output.writeLine(join((contactMap.data[contactMap.width * i .. contactMap.width * (i + 1) - 1]), " "))
output.close()

let endsfile = open("jurkat_binends_100K.txt", fmWrite)
endsfile.writeline("chr\tend")
for key in bins.keys:
    for pos in bins[key]:
        endsfile.writeLine(key, "\t", pos)
endsfile.close()

echo "Done!"
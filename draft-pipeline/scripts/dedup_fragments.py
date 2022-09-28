import sys

# Merge duplicate fragments sorted by (chr, start, end, cell) 
# Fragment should be format chr, start, end, cell, count (where count can be 1)
# In output, the counts of adjacent of (chr, start, end, cell) will be added together
# so only one line is output per unique (chr, start, end, cell) combination

# Authors: Ben Parks
# Last updated: 9/27/22

# Usage: cat input.tsv | python dedup_fragments.py > output.tsv

prev_frag = None
prev_count = 0

# Compare each line to the previous, and add up the duplicate counts for identical lines
for line in sys.stdin.buffer:
    idx = line.rfind(b"\t")
    frag = line[:idx]
    count = int(line[idx+1:])
    if frag == prev_frag:
        prev_count += count
        continue
    if prev_frag:
        sys.stdout.buffer.write(prev_frag)
        sys.stdout.buffer.write(b"\t")
        sys.stdout.buffer.write(str(count).encode())
        sys.stdout.buffer.write(b"\n")
    prev_frag = frag
    prev_count = count
        

# Write out the last line
if prev_frag:
    sys.stdout.buffer.write(prev_frag)
    sys.stdout.buffer.write(b"\t")
    sys.stdout.buffer.write(str(count).encode())
    sys.stdout.buffer.write(b"\n")
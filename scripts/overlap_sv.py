#!/usr/bin/env python


if __name__ == "__main__":
    sv_truth = dict()
    # SV bed file
    with open("simulated_20.bed") as infile:
        for line in infile:
            con, start, con2, stop, mutype = line.strip().split()
            sv_truth[start] = [abs(int(stop) - int(start)), mutype]

    manta_miss = []
    # annotate with missing SV:
    with open("manta_res.tsv") as infile:
        for line in infile:
            if line.startswith("0") or line.startswith("4"):
                manta_miss.append(line.strip().split()[2])

    with open("manta_join.tsv", "w") as outfile:
        for key, value in sv_truth.items():
            if key in manta_miss:
                outfile.write(f"{key} {value[0]} {value[1]}  missing\n")
            else:
                outfile.write(f"{key} {value[0]} {value[1]}  found\n")

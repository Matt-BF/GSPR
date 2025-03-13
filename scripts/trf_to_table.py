import sys

with open(sys.argv[1]) as f:
    result_dict = {}
    line = f.readline()
    while line:
        if line.startswith("@"):
            print(f"working on {line}")
            contig_line = line.strip().replace("@", "")
            result_dict[contig_line] = []
            line = f.readline()
            try:
                while line.split()[0].isdigit():
                    # results_stats = line.split()[0:13]
                    # results_repeats = line.split()[13:]
                    result_dict[contig_line].append(line.split())
                    line = f.readline()
            except IndexError:
                break


with open("test.tsv", "w") as out_table:
    out_table.write(
        "contig_name\trepeat_start\trepeat_end\tperiod_size\tnum_copies\tpattern_size\tpercent_matches\tpercent_indels\talignment_score\ta_composition\tt_composition\tc_composition\tg_composition\tentropy\trepeat\tpattern\n"
    )
    for contig in result_dict:
        for result in result_dict[contig]:
            result_stats = "\t".join(result[0:13])
            result_pattern = result[13:]
            out_table.write(
                f"{contig}\t{result_stats}\t{result_pattern[0]}\t{','.join(result_pattern[1:])}\n"
            )

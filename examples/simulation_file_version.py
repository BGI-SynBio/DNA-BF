#!/usr/bin/env python3
#-*- coding:utf-8 -*-

from bean import bloom
from exps.pipelines import ErrorCounterOneFile
import random


def dna_file(dna_file_path, extract_count=None, seed=None, extracted_file=None):
    digital_dnas = []
    with open(dna_file_path, 'r') as file:
        lines = file.readlines()
        if extract_count:
            if seed:
                random.seed(seed)
            extract_lines = random.sample(lines, extract_count)
            if extracted_file:
                with open(extracted_file + "-" + str(extract_count) + "-" + str(seed) + ".dna", "w") as f:
                    f.writelines(extract_lines)
            for line in extract_lines:
                digital_dnas.append(line[:-1])
        else:
            for line in lines:
                digital_dnas.append(line[:-1])

    return digital_dnas


if __name__ == '__main__':

    # synthesizing and sequencing simulation of all sequences from different versions,
    # and realization of file version control function
    file_lst = ["version1.dna", "version2.dna", "version3.dna", "version4.dna", "version5.dna"]

    all_version = []
    for file in file_lst:
        with open("../files/version/" + str(file), "r") as f:
            lines = f.readlines()
            for line in lines:
                if line[:-1] not in all_version:
                    all_version.append(line[:-1])

    for idx, file in enumerate(file_lst):
        if idx==0:
            one_version = dna_file("../files/version/" + str(file), seed=30)
            training_name = "version-" + str(idx)
            checking_name = "all_version"

            bloom_filter = bloom.Filter(inform_count=len(one_version), hash_function_type="FNV1a")
            for depth in [100]:
                error_counter = ErrorCounterOneFile(bloom_filter, expected_rate=0.001, payload_length=120,
                                                    digital_count=len(one_version), avg_yield=1000, stdev=100,
                                                    read_depth=100, root_path="../outputs/", insertion_rate=0.00075,
                                                    deletion_rate=0.00075, substitution_rate=0.0015, seed=30,
                                                    actual_dnas=True, training_dnas=one_version,
                                                    checking_dnas=all_version, training_name=training_name,
                                                    checking_name=checking_name, filter_size=0, hash_size=0,
                                                    counting_type=True, record_coverage=False, coverage_threshold=None,
                                                    need_record=False, verbose=False)
                error_counter.run()
                error_counter.show_results()

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

    # synthesizing and sequencing simulation of purely random generated sequences,
    # and realization of the anti-contamination function
    bloom_filter = bloom.Filter(inform_count=100000, hash_function_type="FNV1a")
    for depth in [1000]:
        error_counter = ErrorCounterOneFile(bloom_filter, expected_rate=0.001, payload_length=120, digital_count=100000,
                                            avg_yield=1000, stdev=100, read_depth=depth, root_path="../outputs/",
                                            insertion_rate=0.00075, deletion_rate=0.00075, substitution_rate=0.0015,
                                            seed=30, actual_dnas=True, training_dnas=None, checking_dnas=None,
                                            training_name=None, checking_name=None, filter_size=0, hash_size=0,
                                            counting_type=False, record_coverage=False, coverage_threshold=12,
                                            need_record=False, verbose=False)
        error_counter.run()
        error_counter.show_results()

    for depth in list(range(10, 1000, 10)):
        error_counter = ErrorCounterOneFile(bloom_filter, expected_rate=0.001, payload_length=120, digital_count=100000,
                                            avg_yield=1000, stdev=100, read_depth=depth, root_path="../outputs/",
                                            insertion_rate=0.00075, deletion_rate=0.00075, substitution_rate=0.0015,
                                            seed=30, actual_dnas=True, training_dnas=None, checking_dnas=None,
                                            training_name=None, checking_name=None, filter_size=0, hash_size=0,
                                            counting_type=False, record_coverage=False, coverage_threshold=12,
                                            need_record=False, verbose=False)
        error_counter.run()
        error_counter.show_results()

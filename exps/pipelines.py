import copy

from bean import bloom
from exps.operations import synthesising, occur_errors, sequencing
from utils.recorder import Recoder
from utils.recorder_seq import Recoder_Seq
import os
from utils.model_saver import save_model, load_model
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import random
from collections import Counter
import sys


class ErrorCounterOneFile(object):

    def __init__(self, bloom_filter, expected_rate, payload_length, digital_count, avg_yield, stdev, read_depth,
                 root_path, insertion_rate=0.00075, deletion_rate=0.00075, substitution_rate=0.0015, seed=None,
                 actual_dnas=False, training_dnas=None, checking_dnas=None, training_name=None, checking_name=None,
                 filter_size=0, hash_size=0, counting_type=False, record_coverage=False, coverage_threshold=None,
                 need_record=False, verbose=False):
        """
        initialize the error counter.

        :param bloom_filter: un-trained Bloom Filter.
        :type: list

        :param expected_rate: expected rate of the Bloom Filter.
        :type: float

        :param payload_length: payload length in each DNA sequence.
        :type: int

        :param digital_count: digital count of DNA sequences.
        :type: int

        :param avg_yield: average value of yield of each DNA sequence.
        :type: int

        :param stdev: stdev of DNA sequence after DNA synthesising.
        :type: int

        :param read_depth: read depth in DNA sequencing.
        :type: int

        :param insertion_rate: insertion rate.
        :type: float

        :param deletion_rate: deletion rate.
        :type: float

        :param substitution_rate: substitution rate.
        :type: float

        :param root_path: the root path of the files.
        :type: string

        :param seed: random seed.
        :type: int

        :param actual_dnas: whether to use the actual DNA sequences.
        :type: bool

        :param training_dnas: actual DNA sequences for training the Bloom Filter.
        :type: list

        :param checking_dnas: actual DNA sequences for checking the Bloom Filter.
        :type: list

        :param training_name: the name of the training_dnas.
        :type: string

        :param checking_name: the name of the checking_dnas.
        :type: string

        :param filter_size: size of filter list.
        :type: int

        :param hash_size: size of hash functions.
        :type: int

        :param counting_type: whether to use Counting Bloom Filter.
        :type: bool

        :param record_coverage: whether to record the information of the checked DNA sequences checked in each coverage.
        :type: bool

        :param coverage_threshold: the coverage threshold for removing false positive sequences.
        :type: int

        :param need_record: record sequencing data and purely random sequences.
        :type bool

        :param verbose: need to show process.
        :type: bool
        """

        self.avg_yield = avg_yield
        self.payload_length = payload_length
        self.stdev = stdev
        self.read_depth = read_depth
        self.insertion_rate = insertion_rate
        self.deletion_rate = deletion_rate
        self.substitution_rate = substitution_rate

        self.root_path = root_path
        self.seed = seed
        self.actual_dnas = actual_dnas

        self.training_dnas = training_dnas
        self.checking_dnas = checking_dnas
        self.training_name = training_name
        self.checking_name = checking_name

        self.counting_type = counting_type
        self.record_coverage = record_coverage
        self.coverage_threshold = coverage_threshold
        self.verbose = verbose
        self.need_record = need_record

        if not self.actual_dnas:
            self.digital_count = digital_count
            self.recorder = Recoder()
        else:
            self.digital_count = len(training_dnas)
            self.recorder = Recoder(digital_dnas=training_dnas)
        self.checked_keys = {}
        self.deleted_list = None

        self.recorder_seq = Recoder_Seq()

        self.filter = bloom_filter
        if filter_size and hash_size:
            self.filter.create_by_customized(filter_size, hash_size)
        else:
            self.filter.create_by_expected(expected_rate)

        self.decoded_count = 0
        self.error_count = 0
        self.unexpected_count = 0
        self.read_count = 0

        self.limit_decoded_count = 0

        self.max_right = None
        self.min_right = None
        self.max_wrong = None
        self.min_wrong = None

        self.right_check = {}
        self.wrong_check = {}
        self.unexpected_dna = {}

        self.false_positive_lst = []
        self.fp_lst_in_wrong = []
        self.fp_lst_in_right = []

        self.file = None
        self.run_time = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.expected_rate = expected_rate

        self.max_other_right_dnas = None
        self.min_other_right_dnas = None
        self.max_actual_error_dnas = None
        self.min_actual_error_dnas = None

    def run(self):
        """
        run the pipeline.
        """
        if not self.actual_dnas:

            file_marker = "random-generate_" + str(self.digital_count) + "_sd" + str(self.seed) \
                          + "_cp" + str(self.avg_yield) + "-std" + str(self.stdev) \
                          + "_ins" + str(self.insertion_rate) \
                          + "-del" + str(self.deletion_rate) \
                          + "-sub" + str(self.substitution_rate)

        else:

            file_marker = str(self.training_name) + "_sd" + str(self.seed) \
                          + "_cp" + str(self.avg_yield) + "-std" + str(self.stdev)\
                          + "_ins" + str(self.insertion_rate)\
                          + "-del" + str(self.deletion_rate)\
                          + "-sub" + str(self.substitution_rate)

        self.file = self.root_path + file_marker
        pkl_path = self.file + ".pkl"

        T1 = datetime.datetime.now()

        if not os.path.exists(pkl_path):
            # digital synthesis
            if not self.actual_dnas:
                t1 = datetime.datetime.now()
                self.recorder.generated_by_random(self.seed, self.payload_length, self.digital_count)
                t2 = datetime.datetime.now()
                self.run_time[0] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000

                t1 = datetime.datetime.now()
                self.filter.train(self.recorder.get_digital_dnas(), counting_type=self.counting_type,
                                  verbose=self.verbose)
                t2 = datetime.datetime.now()
                self.run_time[1] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000
            else:
                t1 = datetime.datetime.now()
                self.filter.train(self.training_dnas, counting_type=self.counting_type, verbose=self.verbose)
                t2 = datetime.datetime.now()
                self.run_time[1] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000

            # synthesis
            t1 = datetime.datetime.now()
            if self.checking_dnas:
                temp_yields = synthesising(self.checking_dnas, self.avg_yield, self.stdev, self.seed, self.verbose)  # suitable to all version
            else:
                temp_yields = synthesising(self.recorder.get_digital_dnas(), self.avg_yield, self.stdev, self.seed, self.verbose)
            t2 = datetime.datetime.now()
            self.run_time[2] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000

            self.recorder.save_current_pool(temp_yields)

            # occur error
            t1 = datetime.datetime.now()
            temp_yields = occur_errors(self.recorder.get_dna_pool(), self.insertion_rate, self.deletion_rate,
                                       self.substitution_rate, self.seed, self.verbose)
            t2 = datetime.datetime.now()
            self.run_time[3] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000
            self.recorder.save_current_pool(temp_yields)
            save_model(pkl_path, self.recorder)

        else:
            t1 = datetime.datetime.now()
            self.recorder = load_model(pkl_path)
            t2 = datetime.datetime.now()
            self.run_time[0] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000
            print("Load pkl file successfully!", flush=True)

            if not self.actual_dnas:
                t1 = datetime.datetime.now()
                self.filter.train(self.recorder.get_digital_dnas(), counting_type=self.counting_type,
                                  verbose=self.verbose)
                # self.filter.train(self.training_dnas, counting_type=self.counting_type, verbose=self.verbose)
                t2 = datetime.datetime.now()
                self.run_time[1] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000
            else:
                t1 = datetime.datetime.now()
                # self.filter.train(self.recorder.get_digital_dnas(), counting_type=self.counting_type,
                #                   verbose=self.verbose)
                self.filter.train(self.training_dnas, counting_type=self.counting_type, verbose=self.verbose)  # suitable to all version
                t2 = datetime.datetime.now()
                self.run_time[1] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000

        # get sequencing data
        if self.need_record:

            new_sequencing_dnas = sequencing(self.recorder.get_dna_pool(), self.read_depth, self.digital_count,
                                             start_depth=0, seed=self.seed, verbose=self.verbose)

            if not self.actual_dnas:
                with open(self.root_path + "random-generate_" + str(self.digital_count) + "_sd" + str(self.seed)
                          + ".dna", "w") as f:
                    f.writelines([x + "\n" for x in self.recorder.get_digital_dnas()])
                pd.DataFrame(new_sequencing_dnas.items()).to_csv(self.root_path + "random-generate_"
                                                                 + str(self.digital_count) + "_"
                                                                 + "sequencing-data_depth" + str(self.read_depth)
                                                                 + ".csv")
            else:
                df = pd.DataFrame(new_sequencing_dnas.items())
                df.rename(columns={0: "DNA", 1: "count"}, inplace=True)
                df.to_csv(self.root_path + str(self.training_name) + "_" + "sequencing-data_depth"
                          + str(self.read_depth) + ".csv", index=False)
            sys.exit(0)

        #  pickle.dump、pickle.load
        t1 = datetime.datetime.now()
        seq_interval = int(10)

        prefix = '_within'
        suffix = 'depth'
        save_depths = [0]
        all_files = os.listdir(self.root_path)
        pkl_files = [file for file in all_files if file.endswith('.pkl')]
        for p_file in pkl_files:
            if p_file.endswith(suffix + '.pkl'):
                if file_marker + "_rightcheck" + str(self.filter._hash_function_type[0]) + "_fs" \
                        + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix in p_file:
                    save_depths.append(int(p_file.split(suffix + ".pkl")[-2].split(prefix)[-1]))

        save_depths.sort()
        max_depth = max(save_depths)
        if self.read_depth <= max_depth:
            for idx, depth in enumerate(list(range(seq_interval, self.read_depth + 1, seq_interval))):

                current_right_check = load_model(self.file + "_rightcheck" + str(self.filter._hash_function_type[0])
                                                 + "_fs" + str(self.filter._filter_size) + "_hs"
                                                 + str(self.filter._hash_size) + prefix + str(depth) + suffix + ".pkl")

                current_wrong_check = load_model(self.file + "_wrongcheck" + str(self.filter._hash_function_type[0])
                                                 + "_fs" + str(self.filter._filter_size) + "_hs"
                                                 + str(self.filter._hash_size) + prefix + str(depth) + suffix + ".pkl")

                current_unexpected_dna = load_model(self.file + "_unexpecteddna"
                                                    + str(self.filter._hash_function_type[0]) + "_fs"
                                                    + str(self.filter._filter_size) + "_hs"
                                                    + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                    + ".pkl")

                current_read_count = load_model(self.file + "_readcount" + str(self.filter._hash_function_type[0])
                                                 + "_fs" + str(self.filter._filter_size) + "_hs"
                                                 + str(self.filter._hash_size) + prefix + str(depth) + suffix + ".pkl")

                current_error_count = load_model(self.file + "_errorcount" + str(self.filter._hash_function_type[0])
                                                 + "_fs" + str(self.filter._filter_size) + "_hs"
                                                 + str(self.filter._hash_size) + prefix + str(depth) + suffix + ".pkl")

                current_unexpected_count = load_model(self.file + "_unexpectedcount"
                                                      + str(self.filter._hash_function_type[0]) + "_fs"
                                                      + str(self.filter._filter_size) + "_hs"
                                                      + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                      + ".pkl")
                if self.counting_type:
                    current_checked_keys = load_model(self.file + "_checkedkeys"
                                                      + str(self.filter._hash_function_type[0]) + "_fs"
                                                      + str(self.filter._filter_size) + "_hs"
                                                      + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                      + ".pkl")
                    current_deleted_list = load_model(self.file + "_deletedlist"
                                                      + str(self.filter._hash_function_type[0]) + "_fs"
                                                      + str(self.filter._filter_size) + "_hs"
                                                      + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                      + ".pkl")

                start_time = datetime.datetime.now()
                self.right_check = dict(Counter(self.right_check) + Counter(current_right_check))
                self.wrong_check = dict(Counter(self.wrong_check) + Counter(current_wrong_check))
                self.unexpected_dna = dict(Counter(self.unexpected_dna) + Counter(current_unexpected_dna))
                self.read_count += current_read_count
                self.error_count += current_error_count
                self.unexpected_count += current_unexpected_count
                if self.counting_type:
                    self.checked_keys.update(current_checked_keys)
                    self.deleted_list = current_deleted_list
                print(str(depth) + '___adding variables: ', datetime.datetime.now() - start_time, flush=True)

        else:
            if max_depth > 0:
                for idx, depth in enumerate(list(range(seq_interval, max_depth + 1, seq_interval))):
                    current_right_check = load_model(self.file + "_rightcheck" + str(self.filter._hash_function_type[0])
                                                     + "_fs" + str(self.filter._filter_size) + "_hs"
                                                     + str(self.filter._hash_size) + prefix
                                                     + str(depth) + suffix + ".pkl")

                    current_wrong_check = load_model(self.file + "_wrongcheck" + str(self.filter._hash_function_type[0])
                                                     + "_fs" + str(self.filter._filter_size) + "_hs"
                                                     + str(self.filter._hash_size) + prefix
                                                     + str( depth) + suffix + ".pkl")

                    current_unexpected_dna = load_model(self.file + "_unexpecteddna"
                                                        + str(self.filter._hash_function_type[0])+ "_fs"
                                                        + str(self.filter._filter_size) + "_hs"
                                                        + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                        + ".pkl")

                    current_read_count = load_model(self.file + "_readcount" + str(self.filter._hash_function_type[0])
                                                    + "_fs" + str(self.filter._filter_size) + "_hs"
                                                    + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                    + ".pkl")

                    current_error_count = load_model(self.file + "_errorcount" + str(self.filter._hash_function_type[0])
                                                     + "_fs" + str(self.filter._filter_size) + "_hs"
                                                     + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                     + ".pkl")

                    current_unexpected_count = load_model(self.file + "_unexpectedcount"
                                                          + str(self.filter._hash_function_type[0]) + "_fs"
                                                          + str(self.filter._filter_size) + "_hs"
                                                          + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                          + ".pkl")
                    if self.counting_type:
                        current_checked_keys = load_model(self.file + "_checkedkeys"
                                                          + str(self.filter._hash_function_type[0]) + "_fs"
                                                          + str(self.filter._filter_size) + "_hs"
                                                          + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                          + ".pkl")

                        current_deleted_list = load_model(self.file + "_deletedlist"
                                                          + str(self.filter._hash_function_type[0]) + "_fs"
                                                          + str(self.filter._filter_size) + "_hs"
                                                          + str(self.filter._hash_size) + prefix + str(depth) + suffix
                                                          + ".pkl")

                    self.right_check = dict(Counter(self.right_check) + Counter(current_right_check))
                    self.wrong_check = dict(Counter(self.wrong_check) + Counter(current_wrong_check))
                    self.unexpected_dna = dict(Counter(self.unexpected_dna) + Counter(current_unexpected_dna))
                    self.read_count += current_read_count
                    self.error_count += current_error_count
                    self.unexpected_count += current_unexpected_count
                    if self.counting_type:
                        self.checked_keys.update(current_checked_keys)
                        self.deleted_list = current_deleted_list

            if max_depth != self.read_depth:
                if not self.deleted_list:
                    self.deleted_list = copy.deepcopy(self.filter.deleted_filter_list)
                else:
                    pass
                for depth in list(range(max_depth + seq_interval, self.read_depth+1, seq_interval)):
                    start_time = datetime.datetime.now()
                    right_check = {}
                    wrong_check = {}
                    unexpected_dna = {}
                    read_count = 0
                    error_count = 0
                    unexpected_count = 0
                    # checked_keys = {}

                    new_sequencing_dnas = sequencing(self.recorder.get_dna_pool(), depth, self.digital_count,
                                                     start_depth=depth - seq_interval, seed=self.seed,
                                                     verbose=self.verbose)
                    for dna in new_sequencing_dnas.keys():
                        if not self.counting_type:
                            actual_state = self.filter.check(dna, counting_type=self.counting_type)
                        else:
                            actual_state, self.checked_keys, \
                            self.deleted_list = self.filter.check(dna, counting_type=self.counting_type,
                                                                  checked_dna_keys=self.checked_keys,
                                                                  deleted_list=self.deleted_list)
                        expected_state = self.recorder.search_in_original_dna(dna)
                        dna_count = new_sequencing_dnas[dna]
                        read_count += dna_count

                        if actual_state and expected_state:
                            c = right_check.get(dna, 0)
                            right_check[dna] = c + dna_count

                        if actual_state and (not expected_state):
                            error_count += dna_count
                            c = wrong_check.get(dna, 0)
                            wrong_check[dna] = c + dna_count

                        if not expected_state:
                            unexpected_count += dna_count
                            c = unexpected_dna.get(dna, 0)
                            unexpected_dna[dna] = c + dna_count

                        if self.verbose:
                            print("\r" + str(depth - seq_interval) + "~" + str(depth) + ": "
                                  + "do = " + str(read_count)
                                  + ", total = " + str(seq_interval * self.digital_count)
                                  + " (" + str(seq_interval) + " x " + str(self.digital_count) + ")"
                                  + ", error_count = " + str(self.error_count)
                                  # + ", actual decoded oligo count = " + str(len(saved_dnas))
                                  + ", total oligo count = " + str(self.digital_count) + ".", end="", flush=True)
                    print(" ", flush=True)
                    print(str(depth - seq_interval) + "~" + str(depth) + ": "
                          + str((datetime.datetime.now() - start_time).days * 24 * 60
                                + (datetime.datetime.now() - start_time).seconds / 60
                                + (datetime.datetime.now() - start_time).microseconds / 60000000) + 'min', flush=True)

                    self.right_check = dict(Counter(self.right_check) + Counter(right_check))
                    self.wrong_check = dict(Counter(self.wrong_check) + Counter(wrong_check))
                    self.unexpected_dna = dict(Counter(self.unexpected_dna) + Counter(unexpected_dna))
                    self.read_count += read_count
                    self.error_count += error_count
                    self.unexpected_count += unexpected_count

                    if self.counting_type:
                        # self.checked_keys.update(checked_keys)
                        self.recorder_seq.save_current_depth_dnas(right_check, wrong_check, unexpected_dna,
                                                                  read_count, error_count, unexpected_count,
                                                                  depth, check_keys=self.checked_keys,
                                                                  deleted_list=self.deleted_list)
                    else:
                        self.recorder_seq.save_current_depth_dnas(right_check, wrong_check, unexpected_dna,
                                                                  read_count, error_count, unexpected_count,
                                                                  depth)

                    save_model(self.file + "_rightcheck" + str(self.filter._hash_function_type[0])+ "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", right_check)
                    save_model(self.file + "_wrongcheck" + str(self.filter._hash_function_type[0]) + "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", wrong_check)
                    save_model(self.file + "_unexpecteddna" + str(self.filter._hash_function_type[0]) + "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", unexpected_dna)
                    save_model(self.file + "_readcount" + str(self.filter._hash_function_type[0]) + "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", read_count)
                    save_model(self.file + "_errorcount" + str(self.filter._hash_function_type[0]) + "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", error_count)
                    save_model(self.file + "_unexpectedcount" + str(self.filter._hash_function_type[0]) + "_fs"
                               + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                               + str(depth) + suffix + ".pkl", unexpected_count)
                    if self.counting_type:
                        save_model(self.file + "_checkedkeys" + str(self.filter._hash_function_type[0]) + "_fs"
                                   + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                                   + str(depth) + suffix + ".pkl", self.checked_keys)
                        save_model(self.file + "_deletedlist" + str(self.filter._hash_function_type[0]) + "_fs"
                                   + str(self.filter._filter_size) + "_hs" + str(self.filter._hash_size) + prefix
                                   + str(depth) + suffix + ".pkl", self.deleted_list)

        t2 = datetime.datetime.now()
        self.run_time[4] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000

        t1 = datetime.datetime.now()
        self.decoded_count = len(self.right_check)
        for count in self.right_check.values():
            if count > max(self.wrong_check.values(), default=0):
                self.limit_decoded_count += 1
                
        self.max_right = max(self.right_check.values(), default=0)
        self.min_right = min(self.right_check.values(), default=0)
        self.max_wrong = max(self.wrong_check.values(), default=0)
        self.min_wrong = min(self.wrong_check.values(), default=0)

        # record the coverage of non-target sequences(belong to error sequences or other versions) in each version.
        if self.training_dnas != self.checking_dnas and self.training_dnas != None:
            other_right_dnas = set(self.checking_dnas).intersection(set(self.wrong_check.keys()))
            actual_error_dnas = set(self.wrong_check.keys()) - other_right_dnas
            if len(other_right_dnas) > 0:
                self.max_other_right_dnas = max([self.wrong_check[x] for x in other_right_dnas])
                self.min_other_right_dnas = min([self.wrong_check[x] for x in other_right_dnas])
            if len(actual_error_dnas) > 0:
                self.max_actual_error_dnas = max([self.wrong_check[x] for x in actual_error_dnas])
                self.min_actual_error_dnas = min([self.wrong_check[x] for x in actual_error_dnas])

        if self.record_coverage:
            wrong_coverage = {}
            for cov in list(range(self.min_wrong, self.max_wrong + 1)):
                count = list(self.wrong_check.values()).count(cov) * cov
                type_num = list(self.wrong_check.values()).count(cov)
                wrong_coverage[cov] = [count, type_num]

            right_coverage = {}
            for cov in list(range(self.min_right, self.max_right + 1)):
                count = list(self.right_check.values()).count(cov) * cov
                type_num = list(self.right_check.values()).count(cov)
                right_coverage[cov] = [count, type_num]

            df_wrong = pd.DataFrame()
            df_wrong[0] = list(wrong_coverage.keys())
            df_wrong[1] = [x[0] for x in wrong_coverage.values()]
            df_wrong[2] = [x[1] for x in wrong_coverage.values()]
            df_wrong.to_csv(self.root_path + "wrong-coverage_" + file_marker + "_"
                            + str(self.filter._hash_function_type[0]) + "_fs" + str(self.filter._filter_size) + "_hs"
                            + str(self.filter._hash_size) + "_rd" + str(self.read_depth) + ".csv",
                            header=['coverage', 'total_count', 'type_number'], index=None)
            del df_wrong, wrong_coverage

            df_right = pd.DataFrame()
            df_right[0] = list(right_coverage.keys())
            df_right[1] = [x[0] for x in right_coverage.values()]
            df_right[2] = [x[1] for x in right_coverage.values()]
            df_right.to_csv(self.root_path + "right-coverage_" + file_marker
                            + "_" + str(self.filter._hash_function_type[0]) + "_fs" + str(self.filter._filter_size)
                            + "_hs" + str(self.filter._hash_size) + "_rd" + str(self.read_depth) + ".csv",
                            header=['coverage', 'total_count', 'type_number'], index=None)
            del df_right, right_coverage

        t2 = datetime.datetime.now()
        self.run_time[5] = (t2 - t1).days*24*60 + (t2 - t1).seconds/60 + (t2 - t1).microseconds/60000000
        print('Time for making statistics: ', self.run_time[5], flush=True)

        if self.counting_type:
            t1 = datetime.datetime.now()
            trained_dnas = self.recorder.get_digital_dnas() if not self.actual_dnas else self.training_dnas
            self.false_positive_lst = self.filter.find_false_positive(self.checked_keys, deleted_list=self.deleted_list,
                                                                      verbose=self.verbose)
            self.fp_lst_in_wrong = set(self.false_positive_lst).intersection(set(self.wrong_check.keys()))
            self.fp_lst_in_right = set(self.false_positive_lst).intersection(set(self.right_check.keys()))

            t2 = datetime.datetime.now()
            self.run_time[6] = (t2 - t1).days * 24 * 60 + (t2 - t1).seconds / 60 + (t2 - t1).microseconds / 60000000
            # print('Time for finding false positive sequences：', self.run_time[6], flush=True)

        else:
            # removing false positive sequences by coverage threshold
            t1 = datetime.datetime.now()
            if self.coverage_threshold:
                final_dnas = pd.DataFrame([(k, v) for k, v in self.wrong_check.items() if v > self.coverage_threshold]
                                          + [(k, v) for k, v in self.right_check.items() if
                                             v > self.coverage_threshold])
                final_dnas.to_csv(self.file + "_final_dnas-Anti_Contamination.txt", index=None, header=None, sep=",")
                t2 = datetime.datetime.now()
                self.run_time[7] = (t2 - t1).days * 24 * 60 + (t2 - t1).seconds / 60 + (t2 - t1).microseconds / 60000000

            else:
                pass

        T2 = datetime.datetime.now()
        self.run_time[8] = (t2 - t1).days*24*60 + (T2 - T1).seconds/60 + (t2 - t1).microseconds/60000000

    def show_results(self):
        """
        show the distribution, and save statistics.
        """
        decoded_rate = str(round(self.decoded_count / self.digital_count * 100, 2)) + "%"
        false_positive_rate_based_on_count = str(round(self.error_count / self.unexpected_count * 100, 6)) + "%"
        false_positive_rate_based_on_type = str(round(len(self.wrong_check) /
                                                      len(self.unexpected_dna) * 100, 6)) + "%"

        passed_dna = {}
        for key in self.right_check.keys():
            c = passed_dna.get(self.right_check[key], 0)
            passed_dna[self.right_check[key]] = c + 1
        c_dna_sorted = sorted(passed_dna.items(), key=lambda k: k[0])
        x = [i[0] for i in c_dna_sorted]
        y = [i[1] for i in c_dna_sorted]
        passed_dna_count = sum(list(map(lambda a, b: a*b, x, y)))

        with open(self.file + "_" + str(self.filter._hash_function_type[0]) + "_fs" + str(self.filter._filter_size)
                  + "_hs" + str(self.filter._hash_size) + "_rd" + str(self.read_depth) + ".csv", 'w+') as file:
            file.write("filter size" + "," + str(self.filter._filter_size) + "\n"
                       + "hash size" + ',' + str(self.filter._hash_size) + "\n"

                       + "passed dna count" + "," + str(passed_dna_count) + "\n"
                       + "decoded count" + "," + str(self.decoded_count) + "\n"
                       + "limit decoded count" + "," + str(self.limit_decoded_count) + "\n"
                       + "digital count" + "," + str(self.digital_count) + "\n"
                       + "error count" + "," + str(self.error_count) + "\n"
                       + "unexpected count" + "," + str(self.unexpected_count) + "\n"
                       + "unexpected type" + "," + str(len(self.unexpected_dna)) + "\n"

                       + "decoded rate" + "," + str(decoded_rate) + "\n"
                       + "false positive rate (based on count)" + "," + str(false_positive_rate_based_on_count) + "\n"
                       + "false positive rate (based on type)" + "," + str(false_positive_rate_based_on_type) + "\n"
                       + "maximum count in right" + "," + str(self.max_right) + "\n"
                       + "minimum count in right" + "," + str(self.min_right) + "\n"
                       + "maximum count in wrong" + "," + str(self.max_wrong) + "\n"
                       + "minimum count in wrong" + "," + str(self.min_wrong) + "\n"

                       + "max_other_right_dnas" + "," + str(self.max_other_right_dnas) + "\n"
                       + "min_other_right_dnas" + "," + str(self.min_other_right_dnas) + "\n"
                       + "min_other_right_dnas" + "," + str(self.max_actual_error_dnas) + "\n"
                       + "min_actual_error_dnas" + "," + str(self.min_actual_error_dnas) + "\n"

                       + "error type" + "," + str(len(self.wrong_check.keys())) + "\n"

                       + "false positive in wrong(based on type)" + "," + str(len(self.fp_lst_in_wrong)) + "\n"
                       + "false positive in right(based on type)" + "," + str(len(self.fp_lst_in_right)) + "\n"
                       + "total number of find_false_positive " + "," + str(len(set(self.false_positive_lst))) + "\n"

                       + "time for generating DNAs or load pkl files" + "," + str(self.run_time[0]) + "min\n"
                       + "time for training BloomFilter" + "," + str(self.run_time[1]) + "min\n"
                       + "time for synthesizing" + "," + str(self.run_time[2]) + "min\n"
                       + "time for occurring error" + "," + str(self.run_time[3]) + "min\n"
                       + "time for sequencing and check BloomFilter" + "," + str(self.run_time[4]) + "min\n"

                       + "time for making statistics" + "," + str(self.run_time[5]) + "min\n"
                       + "time for finding false positive sequences" + "," + str(self.run_time[6]) + "min\n"
                       + "time for removing false positive sequences by coverage threshold" + ","
                       + str(self.run_time[7]) + "min\n"

                       + "total time" + "," + str(self.run_time[8]) + "min")

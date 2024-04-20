import math
import bean.hasher as hasher
import numpy as np
import copy


class Filter(object):

    def __init__(self, inform_count, hash_function_type):
        """
        Initialize the Bloom Filter through the hash size and bit array size.

        :param inform_count: count of inform.
        :type: int

        :param hash_function_type: type of hash function.
        :type: string
        """
        if inform_count <= 0:
            raise ValueError("wrong data count: " + str(inform_count) + ".")

        self._data_count = inform_count
        self._hash_function_type = hash_function_type

        self._filter_list = None
        self._filter_size = 0
        self._hash_size = 0
        self._is_trained = False

        self.deleted_filter_list = None

    def create_by_expected(self, expected_rate):
        """
        Create the specific hyper-parameters by expected false positive rate.

        :param expected_rate: expected false positive rate.
        :type: float
        """
        if expected_rate <= 0 or expected_rate >= 1:
            raise Exception("Input expected false positive rate is: " + str(expected_rate) +
                            ", expected false positive rate should larger than 0 and less than 1.")

        self._filter_size = int(math.ceil(-self._data_count * math.log(expected_rate) / math.pow(math.log(2), 2)))
        # print("filter size = " + str(self._filter_size), flush=True)
        if self._filter_size % 8 != 0:
            self._filter_size += 8 - self._filter_size % 8
        self._filter_list = [0 for _ in range(self._filter_size)]
        self._hash_size = int(math.ceil(self._filter_size / self._data_count * math.log(2)))

    def create_by_customized(self, filter_size=0, hash_size=0):
        """
        Create the specific hyper-parameters by expected false positive rate.

        :param filter_size: size of filter list.
        :type: int

        :param hash_size: size of hash functions.
        :type: int
        """
        if filter_size <= 0:
            raise ValueError("wrong filter size: " + str(filter_size) + ".")

        self._filter_size = filter_size
        if self._filter_size % 8 != 0:
            self._filter_size += 8 - self._filter_size % 8
        self._filter_list = [0 for _ in range(self._filter_size)]

        if hash_size <= 0:
            self._hash_size = int(math.ceil(self._filter_size / self._data_count * math.log(2)))
        else:
            self._hash_size = hash_size

    def train(self, fragments, counting_type=False, saved_arr=None, verbose=False):
        """
        Train the Bloom Filter by informs.

        :param fragments: all inform for filter.
        :type fragments: list

        :param counting_type: whether to use Counting Bloom Filter.
        :type: bool

        :param saved_arr: the pathname of the saved BF or CBF.
        :type: string or none

        :param verbose: need to show process.
        :type verbose: bool
        """
        if len(fragments) > self._data_count:
            raise ValueError("the self._data_count " + str(self._data_count) +
                             " must be grater than count of inputted DNA sequences " + str(len(fragments)) + ".")

        if saved_arr:
            self._filter_list = list(np.load(saved_arr))
            self._is_trained = True
            if counting_type:
                self.deleted_filter_list = copy.deepcopy(self._filter_list)

        else:
            if verbose:
                print("train Bloom Filter:", flush=True)
            fragments = set(fragments)
            for index, fragment in enumerate(fragments):
                if verbose:
                    print("\rdo = " + str(index + 1) + ", total = " + str(len(fragments)) + ".", end="", flush=True)

                for hash_index in range(self._hash_size):
                    # key = hasher.function(self._hash_function_type[hash_index % len(self._hash_function_type)],
                    #                       fragment, hash_index) % self._filter_size
                    key = hasher.function(self._hash_function_type, fragment, hash_index) % self._filter_size
                    if not counting_type:
                        self._filter_list[key] = 1
                    else:
                        self._filter_list[key] += 1

            if counting_type:
                self.deleted_filter_list = copy.deepcopy(self._filter_list)
                print("Train __ original length of self.deleted_filter_list",
                      len(np.where(np.array(self.deleted_filter_list) < 0)[0]))
                print("length of self.deleted_filter_list", len(self.deleted_filter_list))
                print("(self._filter_list) > 0", len(np.where(np.array(self._filter_list) > 0)[0]))
                print(self.deleted_filter_list)
                print(sum(self.deleted_filter_list))

            print('', flush=True)
            print('max(self._filter_list)', max(self._filter_list), flush=True)
            # if verbose:
            #    print("training finished.")
            print("training finished.", flush=True)

            self._is_trained = True

    def check(self, fragment, counting_type=False, checked_dna_keys=None, deleted_list=None):
        """
        Check one inform by trained filter.

        :param fragment: one requested information fragment.
        :type: string

        :param counting_type: whether to use Counting Bloom Filter.
        :type: bool

        :param checked_dna_keys: record the keys of all the checked fragments.
        :type: dict

        :param deleted_list: bit array has performed deletion operation.
        :type: list

        :return: whether the inform is contained in the trained filter.
        :type: bool
        """
        if not self._is_trained:
            raise Exception("cannot work if the filter is not trained!")

        if not counting_type:
            for hash_index in range(self._hash_size):
                # key = hasher.function(self._hash_function_type[hash_index % len(self._hash_function_type)], fragment,
                #                       hash_index) % self._filter_size
                key = hasher.function(self._hash_function_type, fragment, hash_index) % self._filter_size
                if self._filter_list[key] == 0:
                    return False
            return True

        else:
            # self.deleted_filter_list = copy.deepcopy(self._filter_list)
            keys_lst = []
            state = True
            for hash_index in range(self._hash_size):
                # key = hasher.function(self._hash_function_type[hash_index % len(self._hash_function_type)], fragment,
                #                       hash_index) % self._filter_size
                key = hasher.function(self._hash_function_type, fragment, hash_index) % self._filter_size
                keys_lst.append(key)
                if self._filter_list[key] <= 0:
                    state = False

            if deleted_list:
                if state:
                    if fragment not in checked_dna_keys.keys():
                        checked_dna_keys[fragment] = keys_lst
                        for k in keys_lst:
                            deleted_list[k] -= 1
            else:
                if state:
                    if fragment not in checked_dna_keys.keys():
                        checked_dna_keys[fragment] = keys_lst
                        for k in keys_lst:
                            self.deleted_filter_list[k] -= 1

            return state, checked_dna_keys, deleted_list

    def find_false_positive(self, checked_dna_keys, deleted_list=None, verbose=False):
        """
        Find out all the false positive DNA fragments from the Counting Bloom Filter.

        :param checked_dna_keys: record the keys of all the checked fragments.
        :type: dict

        :param deleted_list: bit array has performed deletion operation.
        :type: list

        :param verbose: need to show process.
        :type verbose: bool

        :return: all the possible false positive dna fragments.
        :rtype: list
        """

        if deleted_list:
            false_positives_keys = np.where(np.array(deleted_list) < 0)[0]
        else:
            false_positives_keys = np.where(np.array(self.deleted_filter_list) < 0)[0]

        if verbose:
            print("Find out all dna fragments whose whole hash values locate at the negative vector.", flush=True)

        false_positive_dnas = []
        fp_keys_set = set(false_positives_keys)
        for index, (dna, keys_lst) in enumerate(checked_dna_keys.items()):
            keys_set = set(keys_lst)
            if fp_keys_set.intersection(keys_set) == keys_set:
                false_positive_dnas.append(dna)
            if verbose:
                print("\rdo = " + str(index + 1) + ", total = " + str(len(checked_dna_keys)), end=" ", flush=True)

        return false_positive_dnas

    def __str__(self):
        return 'BLOOM-FILTER::HFT:%s;HS:%d;FS:%d;BF:%s' % \
               (
                   self._hash_function_type,
                   self._hash_size,
                   self._filter_size,
                   str(self._filter_list)
               )

import copy
import datetime
import random
import numpy as np
from collections import Counter


def synthesising(dnas, avg_yield, stdev, seed=None, verbose=False):
    """
    Synthesising yield of DNA sequences.

    :param dnas: generated DNA sequences.
    :type: list

    :param avg_yield: average value of yield of each DNA sequence.
    :type: int

    :param stdev: standard deviation of the number of each DNA sequence.
    :type: int

    :param seed: random seed.
    :type: int

    :param verbose: need to show process.
    :type: bool

    :return: DNA sequences with yield.
    :rtype: dict
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    print(" ", flush=True)
    print("write DNA sequences: ", flush=True)

    counts = list(map(int, np.random.normal(loc=avg_yield, scale=stdev, size=len(dnas))))
    remaining = avg_yield * len(dnas) - sum(counts)
    for _ in range(abs(remaining)):
        counts[random.randint(0, len(dnas) - 1)] += (1 if remaining > 0 else -1)
        
    index = 0
    dna_yields = {}
    for dna, count in zip(dnas, counts):
        dna_yields[dna] = count

        if verbose:
            print("\rdo = " + str(index + 1) + ", total = " + str(len(dnas)) + ".", end=" ", flush=True)
            index += 1

    print('', flush=True)

    return dna_yields


def sequencing(dna_yields, read_depth, digital_count, start_depth=0, seed=None, verbose=False):
    """
    Sequencing DNA sequence one by one.

    :param dna_yields: DNA sequences with yield.
    :type: dict

    :param read_depth: maximum sequencing count.
    :type: int

    :param digital_count: oligo count generated from file.
    :type: int

    :param start_depth: the initial read depth.
    :type: int

    :param seed: random seed.
    :type: int

    :param verbose: need to show process.
    :type: bool

    :return: current sequenced DNA sequence and its yield.
    """

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    print(" ", flush=True)
    print("Integrate index counts of dna_yields", flush=True)
    dnas = list(dna_yields.keys())
    pool = []
    for index, (dna, count) in enumerate(dna_yields.items()):
        # dnas.append(dna)
        pool += [index] * count

        if verbose:
            print("\rdo = " + str(index + 1) + ", total = " + str(len(dna_yields)) + ".", end=" ", flush=True)

    st = datetime.datetime.now()
    # random.shuffle(pool)
    pool = np.random.permutation(pool)    ## faster
    print(" ", flush=True)
    print("Time for disrupting the order of the pool：", datetime.datetime.now() - st, flush=True)
    print(datetime.datetime.now(), flush=True)

    print(" ", flush=True)
    print("read DNA sequences: ", flush=True)
    # for index in range(read_depth * digital_count):
    #     if verbose:
    #         print("\rdo = " + str(index + 1) + ", total = " + str(read_depth * digital_count) + ".", end=" ", flush=True)
    #
    #     if index < len(pool):
    #         yield dnas[pool[index]]
    #     else:
    #         yield dnas[random.sample(pool, 1)[0]]
    #         # yield dnas[random.sample(pool.tolist(), 1)[0]]

    total_count = read_depth * digital_count
    start_count = start_depth * digital_count
    pool_index = dict(Counter(pool[start_count:total_count]))

    # return all the dnas.
    sequencing_dna_yields = {}
    for idx, index in enumerate(pool_index.keys()):
        # sequencing_dna_yields[list(dna_yields.keys())[index]] = pool_index[index]
        sequencing_dna_yields[dnas[index]] = pool_index[index]

        if verbose:
            print("\rdo = " + str(idx + 1) + ", total  = " + str(len(pool_index)) + ".", end=" ", flush=True)

    return sequencing_dna_yields

    # # yield one dna with its count
    # for idx, index in enumerate(pool_index.keys()):
    #     # dna = list(dna_yields.keys())[index]
    #     dna = dnas[index]
    #     dna_count = pool_index[index]
    #
    #     if verbose:
    #         print("\rdo = " + str(idx + 1) + ", total  = " + str(len(pool_index)) + ".", end=" ", flush=True)
    #
    #     yield dna, dna_count


def occur_errors(dna_yields, insert_rate, delete_rate, replace_rate, seed=None, verbose=False):
    """
    Occur errors in DNA sequences.

    :param dna_yields: DNA sequences with yield.
    :type: dict

    :param insert_rate: the insertion rate  of synthesized DNA.
    :type: float

    :param delete_rate: the deletion rate  of synthesized DNA.
    :type: float

    :param replace_rate: the replacement rate  of synthesized DNA.
    :type: float

    :param seed: random seed.
    :type: int

    :return: current synthesized DNA sequence with unexpected DNA and its yield.
    """
    if seed is not None:
        random.seed(seed)

    dnas = list(dna_yields.keys())
    counts = list(dna_yields.values())
    dnas_set = set(dnas)

    novel_dnas = {}

    print(" ", flush=True)
    print("occur errors:", flush=True)
    sum_modify_count = 0
    for index in range(len(dna_yields)):
        modify_count = 0
        initial_dna = dnas[index]
        for _ in range(counts[index]):
            is_modify = False
            dna = list(copy.deepcopy(initial_dna))

            final_dna = []
            for base in dna:
                final_dna.append(base)
                if random.random() <= insert_rate:
                    is_modify = True
                    final_dna.append(random.choice("ACGT"))
                if random.random() <= replace_rate:
                    is_modify = True
                    final_dna[-1] = random.choice("ACGT".replace(final_dna[-1], ""))
                if random.random() <= delete_rate:
                    is_modify = True
                    del final_dna[-1]

            if is_modify:
                modify_count += 1

                dna = "".join(final_dna)

                if dna in dnas_set:
                    dna_yields[dna] += 1
                elif dna in novel_dnas:
                    novel_dnas[dna] += 1
                else:
                    novel_dnas[dna] = 1

        dna_yields[initial_dna] -= modify_count
        sum_modify_count += modify_count
        if verbose:
            print("\rdo = " + str(index + 1) + ", total = " + str(len(counts)) + ".", end=" ", flush=True)
            print("\rdo = " + str(index + 1) + ", total = " + str(len(dna_yields)) + ",sum of the counts = "
                  + str(sum(list(dna_yields.values()))) + ".", end=" ", flush=True)

    if verbose:
        print("\rsum_modify_count = " + str(sum_modify_count) + ", len(novel_dnas) = " + str(len(novel_dnas))
              + ", sum(novel_counts) = " + str(sum(list(novel_dnas.values()))) + '.\n' + ", len(counts) = "
              + str(len(dna_yields)) + ", sum(counts) = " + str(sum(list(dna_yields.values()))) + ", len(dnas) = "
              + str(len(dna_yields))+'.\n', end=" ", flush=True)

    new_dna_yields = copy.deepcopy(dna_yields)
    new_dna_yields.update(novel_dnas)

    dna_yields = {}

    for dna, count in new_dna_yields.items():
        if count > 0:
            dna_yields[dna] = count

    del new_dna_yields

    if verbose:
        print("sum(counts)_after error occur = ", sum(dna_yields.values()), flush=True)
        print("number of all types:", len(dna_yields), flush=True)

    return dna_yields
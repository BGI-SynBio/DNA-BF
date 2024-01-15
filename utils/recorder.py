import copy
import random
import math


class Recoder(object):

    def __init__(self, digital_dnas=None):
        """
        Initialize the DNA recorder.

        :param digital_dnas: digital DNA sequences encoded by the digital file(s).
        """
        self.digital_dnas = digital_dnas
        self.saved_dna_pool = None
        self.saved_full_depth_index = None
        # self.total_number = 0

    def generated_by_random(self, seed, payload_length, oligo_count, file_count=1):
        """
        Generate original DNA sequences by random.

        :param seed: random seed.
        :param payload_length: length of DNA payload in one file.
        :param oligo_count: oligo count generated from one file.
        :param file_count: number of digital file.
        """
        random.seed(seed)

        self.digital_dnas = []

        index_length = 20
        base_dict = {'00': "A", '01': "C", '10': "G", '11': "T"}
        # index_length = math.ceil(len(bin(oligo_count)[2:]) / 2.0)
        for oligo_index in range(oligo_count):
            number = bin(oligo_index)[2:]
            index_string = ''
            if len(number) % 2 == 0:
                pass
            else:
                number = "0" + number
            for idx in range(0, len(number), 2):
                x = number[idx:idx + 2]
                index_string = index_string + base_dict[x]

            index_string = "A" * (index_length - int(len(number)/2)) + index_string

            dna = [s for s in index_string] + [random.choice("ACGT") for _ in range(payload_length)]
            self.digital_dnas.append("".join(dna))

    def get_digital_dnas(self):
        """
        Obtain digital sequences from a digital file.

        :return: digital sequences from a requested digital file.
        """
        return self.digital_dnas

    def get_dna_pool(self):
        """
        Obtain current DNA pool.

        :return: saved DNA yields.
        """
        return self.saved_dna_pool

    def get_full_depth_index(self):
        """
        Obtain current full depth sequencing index.

        :return: saved full depth sequencing index.
        """
        return self.saved_full_depth_index

    def save_current_pool(self, dna_pool):
        """
        Save current DNA yields

        :param dna_pool: current DNA pool.
        """
        self.saved_dna_pool = dna_pool

    def save_current_full_depth_index(self, index_lst):
        """
        Save current full depth sequencing index.

        :param index_lst: current full depth sequencing index.
        """
        self.saved_full_depth_index = index_lst

    def search_in_original_dna(self, dna):
        """
        Search a DNA sequence from digital DNA sequences.

        :param dna: the searched DNA sequence.

        :return: whether it included in one file.
        """
        return dna in self.digital_dnas

    def search_in_current_dna(self, dna):
        """
        Search a DNA sequence from the saved DNA yields.

        :param dna: the searched DNA sequence.

        :return: the copy number of the searched DNA sequence.
        """
        copy_number = self.saved_dna_pool.get(dna)

        if copy_number is None:
            return 0

        return copy_number

    def get_wrong_rate(self, total_number):
        """
        Obtain wrong rate in the global sequences.

        :param total_number: total number of DNA sequences.

        :return: wrong rate.
        """
        wrong_number = 0

        for dna, count in self.saved_dna_pool.items():
            if dna not in self.digital_dnas:
                wrong_number += count

        return wrong_number / total_number

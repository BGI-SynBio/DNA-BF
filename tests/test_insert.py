import random
import unittest
import math
import bean.bloom as bloom


class TestBloomInsert(unittest.TestCase):

    def setUp(self):
        random.seed(30)
        self.dna_length = 10
        self.count = 10000
        self.false_positive = 0.001
        self.base_list = ["A", "C", "G", "T"]
        self._init_dnas()
        self.filter = bloom.Filter(inform_count=int(math.pow(4, self.dna_length)), hash_function_type=["BKDR"])
        self.filter.create_by_expected(self.false_positive)
        self.filter.train(self.dnas)

    def _init_dnas(self):
        self.dnas = []
        while len(self.dnas) < self.count:
            dna = [random.choice(["A", "C", "G", "T"]) for _ in range(self.dna_length)]
            if dna not in self.dnas:
                self.dnas.append(dna)

    def test_correct_matrix(self):
        results = []
        for fragment in self.dnas:
            results.append(self.filter.check(fragment))

        self.assertEqual(results, [True for _ in range(len(self.dnas))])

    def test_wrong_matrix(self):
        total = int(math.pow(4, self.dna_length + 1))
        error_count = 0
        for number in range(total):
            dna = []
            while number > 0:
                dna.append(self.base_list[number % 4])
                number //= 4
            dna = ["A" for _ in range(self.dna_length + 1 - len(dna))] + dna

            if dna not in self.dnas and self.filter.check(dna):
                error_count += 1

        self.assertLess(error_count / int(math.pow(4, self.dna_length)), self.false_positive)

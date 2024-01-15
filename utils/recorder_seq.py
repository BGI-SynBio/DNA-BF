import copy
import random
import math


class Recoder_Seq(object):

    def __init__(self, depth=0):
        """
        Initialize the recorder for the sequencing DNAs.

        :param depth: the read depth of this experiment.
        :type : int
        """
        self.r = locals()
        self.w = locals()
        self.unexpected = locals()
        self.read_count = locals()
        self.error_count = locals()
        self.unexpected_count = locals()
        self.checked_keys = locals()
        self.deleted_list = locals()
        # for depth in list(range(0+depth_interval, read_depth+1, depth_interval)):
        self.r['self.saved_right_' + str(depth)] = None
        self.w['self.saved_wrong_' + str(depth)] = None
        self.unexpected['self.saved_unexpected_' + str(depth)] = None
        self.read_count['self.saved_read_count_' + str(depth)] = None
        self.error_count['self.saved_error_count_' + str(depth)] = None
        self.unexpected_count['self.saved_unexpected_count_' + str(depth)] = None
        self.checked_keys['self.saved_checked_keys_' + str(depth)] = None
        self.deleted_list['self.saved_deleted_list_' + str(depth)] = None

    def get_one_depth_dnas(self, one_depth):
        """
        Obtain check results under one read depth.

        :param one_depth: current read depth.
        :type: int

        :return: right checked and wrong checked dnas and their counts, unexpected dnas and their counts.
        :rtype: dict
        """
        return self.r['self.saved_right_' + str(one_depth)], self.w['self.saved_wrong_' + str(one_depth)], \
               self.unexpected['self.saved_unexpected_' + str(one_depth)],\
               self.read_count['self.saved_read_count_' + str(one_depth)], \
               self.error_count['self.saved_error_count_' + str(one_depth)], \
               self.unexpected_count['self.saved_unexpected_count_' + str(one_depth)], \
               self.checked_keys['self.saved_checked_keys_' + str(one_depth)], \
               self.deleted_list['self.saved_deleted_list_' + str(one_depth)]


    def save_current_depth_dnas(self, right_check, wrong_check, unexpected_dnas, read_count, error_count,
                                unexpected_count, one_depth, check_keys=None, deleted_list=None):
        """
        Save current depth right checked and wrong checked dnas and their counts.

        :param right_check: right checked dnas and their dna counts.
        :type: dict

        :param wrong_check: wrong checked dnas and their dna counts.
        :type: dict

        :param one_depth: current read depth.
        :type: int
        """
        self.r['self.saved_right_' + str(one_depth)] = copy.deepcopy(right_check)
        self.w['self.saved_wrong_' + str(one_depth)] = copy.deepcopy(wrong_check)
        self.unexpected['self.saved_unexpected_' + str(one_depth)] = copy.deepcopy(unexpected_dnas)

        self.read_count['self.saved_read_count_' + str(one_depth)] = copy.deepcopy(read_count)
        self.error_count['self.saved_error_count_' + str(one_depth)] = copy.deepcopy(error_count)
        self.unexpected_count['self.saved_unexpected_count_' + str(one_depth)] = copy.deepcopy(unexpected_count)

        self.checked_keys['self.saved_checked_keys_' + str(one_depth)] = copy.deepcopy(check_keys)
        self.deleted_list['self.saved_deleted_list_' + str(one_depth)] = copy.deepcopy(deleted_list)

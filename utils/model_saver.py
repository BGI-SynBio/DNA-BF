import sys
import pickle

from utils import log


# noinspection PyProtectedMember
def save_model(path, model, need_log=False):
    """
    introduction: Save model to file.

    :param path: The path of file.
                  Usually in the models directory.

    :param model: Current model for encoding.
                   Type: .pkl
                   e.g. YYC.

    :param need_log: choose to output log file or not.
    """
    if need_log:
        log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                   "Save model to file: " + path)

    with open(path, "wb") as file:
        pickle.dump(model, file)


# noinspection PyProtectedMember
def load_model(path, need_log=False):
    """
    introduction: Load model from file.

    :param path: The path of file.
                  Type: .pkl

    :return: needed model.
              e.g. YYC.

    :param need_log: choose to output log file or not.
    """
    if need_log:
        log.output(log.NORMAL, str(__name__), str(sys._getframe().f_code.co_name),
                   "Load model from file: " + path)

    with open(path, "rb") as file:
        return pickle.load(file)

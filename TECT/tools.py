import os


def mkdir(path):
    folder_exist = os.path.exists(path)
    if not folder_exist:
        os.makedirs(path)

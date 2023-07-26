import os


def load_paths():
    cwd = os.getcwd()
    io_path = os.path.join(cwd, 'IO')

    return io_path

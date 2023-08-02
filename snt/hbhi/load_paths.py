import os


def load_paths(user_path):
    if not user_path:
        user_path = os.path.expanduser('~')

    home_path = user_path
    io_path = os.path.join(home_path, 'IO')

    return io_path

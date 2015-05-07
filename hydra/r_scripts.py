import os.path


def get_script_path(filename):
    '''return full path of R script stored with rnaseq_pipeline package'''
    package_dir = os.path.dirname(__file__)
    return os.path.join(package_dir, 'r_scripts', filename)

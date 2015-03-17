import os.path


def get_script_path(filename):
    '''return full path of R script stored with rnaseq_pipeline package'''
    return os.path.join(__file__, 'r_scripts', filename)

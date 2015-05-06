"""Scripts that handle the logging. Rather simplistic at this stage
"""

def writeLog(string, param):
    """Writes a sring into the log file

    :Parameter string: string that is written to the log file
    :Parameter param: dictionary that contains all general RNASeq pipeline parameters
    """
    if param['verbose']:
        print string
    handle = open(param['log_handle'], 'a')
    handle.write(string)
    handle.close()


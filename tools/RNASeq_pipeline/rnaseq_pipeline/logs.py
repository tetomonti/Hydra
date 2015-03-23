
def writeLog(string,param):
    if param['verbose']: print string
    handle=open(param['log_handle'],'a')
    handle.write(string)
    handle.close()


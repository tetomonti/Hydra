#Copyright 2015 Daniel Gusenleitner, Stefano Monti

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

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


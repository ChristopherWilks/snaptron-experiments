#!/usr/bin/env python2.7

# Copyright 2016, Christopher Wilks <broadsword@gmail.com>
#
# This file is part of Snaptron.
#
# Snaptron is free software: you can redistribute it and/or modify
# it under the terms of the 
# Creative Commons Attribution-NonCommercial 4.0 
# International Public License ("CC BY-NC 4.0").
#
# Snaptron is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# CC BY-NC 4.0 license for more details.
#
# You should have received a copy of the CC BY-NC 4.0 license
# along with Snaptron.  If not, see 
# <https://creativecommons.org/licenses/by-nc/4.0/legalcode>.

import sys
import subprocess
from SnaptronIterator import SnaptronIterator
import clsnapconf

SCRIPT_PATH='../'

class SnaptronIteratorLocal(SnaptronIterator):

    def __init__(self,query_param_string,instance,endpoint):
        SnaptronIterator.__init__(self,query_param_string,instance,endpoint) 
        
        self.construct_query_string(query_param_string,endpoint)
        self.execute_query_string(self.cmd)
    
    def construct_query_string(self,query_param_string,endpoint):
        self.cmd = ['python','%s%s' % (SCRIPT_PATH,clsnapconf.LOCAL_ENDPOINTS[endpoint]), query_param_string]
        return self.cmd

    def execute_query_string(self,query_string):
        #sys.stderr.write("executing %s\n" % (query_string))
        self.subp = subprocess.Popen(query_string,shell=False,stderr=subprocess.PIPE,stdout=subprocess.PIPE,stdin=None)
        self.response = self.subp.stdout
        self.errors = self.subp.stderr
        return self.response
    
    def fill_buffer(self):
        #extend parent version to make sure we close gracefully
        lines_read = SnaptronIterator.fill_buffer(self)
        if lines_read > 0:
            return lines_read
        self.subp.wait()
        errors = self.errors.read()
        #dump any stderr output, not necessarily an error
        #and it if its they'll see the stacktrace from the local module
        if len(errors) > 0:
            sys.stderr.write(errors)
            #raise RuntimeError("error from local command call %s" % (self.cmd))
        return 0

if __name__ == '__main__':
    sIT = SnaptronIteratorLocal('regions=chr2:29446395-30142858&contains=1&rfilter=samples_count>:100&rfilter=annotated:1','snaptron')
    for record in sIT:
        print (record)


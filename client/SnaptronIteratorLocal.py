#!/usr/bin/env python2.7

# This file is part of Snaptron.
#
# Snaptron is free software: you can redistribute it and/or modify
# it under the terms of the 
#
#    The MIT License
#
#    Copyright (c) 2016-  by Christopher Wilks <broadsword@gmail.com> 
#                         and Ben Langmead <langmea@cs.jhu.edu>
#
#    Permission is hereby granted, free of charge, to any person obtaining
#    a copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to
#    permit persons to whom the Software is furnished to do so, subject to
#    the following conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

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


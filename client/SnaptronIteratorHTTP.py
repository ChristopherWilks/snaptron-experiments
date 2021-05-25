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
import urllib2
import httplib
from SnaptronIterator import SnaptronIterator
import clsnapconf
import clsnaputil


class SnaptronIteratorHTTP(SnaptronIterator):

    def __init__(self,query_param_string,instance,endpoint):
        if endpoint == 'sample' or endpoint == 'annotation':
            endpoint += 's'
        SnaptronIterator.__init__(self,query_param_string,instance,endpoint) 

        self.SERVICE_URL=clsnapconf.SERVICE_URL
        self.construct_query_string()
        self.execute_query_string()

    def construct_query_string(self):
        self.query_string = "%s/%s/%s?%s" % (self.SERVICE_URL,self.instance,self.endpoint,self.query_param_string)
        return self.query_string
    
    @clsnaputil.retry((urllib2.HTTPError,urllib2.URLError), tries=17, delay=2, backoff=2)
    def urlopen(self):
       return urllib2.urlopen(self.query_string)

    def execute_query_string(self):
        sys.stderr.write("%s\n" % (self.query_string))
        self.response = self.urlopen()
        return self.response

    def fill_buffer(self):
        #extend parent version to catch HTTP specific error
        try:
            return SnaptronIterator.fill_buffer(self)
        except httplib.IncompleteRead, ir:
            sys.stderr.write(ir.partial)
            raise ir

if __name__ == '__main__':
    it = SnaptronIteratorHTTP('regions=chr1:1-100000&rfilter=samples_count>:5', 'srav1', 'snaptron')
    for r in it:
        print("%s" % r)

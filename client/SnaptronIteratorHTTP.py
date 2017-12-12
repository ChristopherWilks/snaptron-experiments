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
import urllib2
import httplib
from SnaptronIterator import SnaptronIterator
import clsnapconf


class SnaptronIteratorHTTP(SnaptronIterator):

    def __init__(self,query_param_strings,instances,endpoints):
        endpoints_ = []
        for endpoint in endpoints:
            if endpoint == 'sample' or endpoint == 'annotation':
                endpoints_.append(endpoint + 's')
            else:
                endpoints_.append(endpoint)
        endpoints = endpoints_
        SnaptronIterator.__init__(self,query_param_strings,instances,endpoints) 

        self.SERVICE_URL=clsnapconf.SERVICE_URL
        self.construct_query_string()
        self.execute_query_string()

    def construct_query_string(self):
        self.query_string = "%s/%s/%s?%s" % (self.SERVICE_URL,self.instances[self.query_idx],self.endpoints[self.query_idx],self.query_param_strings[self.query_idx])
        return self.query_string

    def execute_query_string(self):
        sys.stderr.write("%s\n" % (self.query_string))
        self.response = urllib2.urlopen(self.query_string)
        return self.response
    
    def next_query_ok(self):
        return self.query_idx < len(self.query_param_strings)

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

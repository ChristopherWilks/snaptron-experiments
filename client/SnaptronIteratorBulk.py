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
import urllib
import urllib2
import httplib
import base64
import time
from SnaptronIterator import SnaptronIterator
import clsnapconf
import clsnaputil


class SnaptronIteratorBulk(SnaptronIterator):

    def __init__(self,query_param_string,instance,endpoint,outfile_handle,processor=None):
        SnaptronIterator.__init__(self,query_param_string,instance,endpoint,processor=processor) 

        self.outfile_handle = outfile_handle
        self.SERVICE_URL=clsnapconf.SERVICE_URL
        self.construct_query_string()
        self.execute_query_string()

    def construct_query_string(self):
        super_string = clsnapconf.BULK_QUERY_DELIMITER.join(self.query_param_string)
        self.data_string = urllib.urlencode({"groups":base64.b64encode(super_string)})
        self.query_string = "%s/%s/%s" % (self.SERVICE_URL,self.instance,clsnapconf.HTTP_ENDPOINTS[self.endpoint])
        return (self.query_string, self.data_string)
    
    @clsnaputil.retry((urllib2.HTTPError,urllib2.URLError), tries=17, delay=2, backoff=2)
    def urlopen(self):
       return urllib2.urlopen(url=self.query_string, data=self.data_string)

    def execute_query_string(self):
        sys.stderr.write("%s\n" % (self.query_string))
        self.response = self.urlopen()
        if self.processor is not None or self.outfile_handle is not None:
            buf_ = self.response.read(self.buffer_size)
            while(buf_ != None and buf_ != ''):
                if self.processor:
                    buf_ += self.response.readline()
                    [self.processor.process(line) for line in buf_.split('\n')]
                else:
                    #since we're doing bulk, just write out as fast as we can
                    self.outfile_handle.write(buf_)
                buf_ = self.response.read(self.buffer_size)

    def fill_buffer(self):
        #extend parent version to catch HTTP specific error
        try:
            return SnaptronIterator.fill_buffer(self)
        except httplib.IncompleteRead, ir:
            sys.stderr.write(ir.partial)
            raise ir

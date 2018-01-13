#!/usr/bin/env python2.7

# Copyright 2017, Christopher Wilks <broadsword@gmail.com>
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
        #since we're doing bulk, just write out as fast as we can
        buf_ = self.response.read(self.buffer_size)
        while(buf_ != None and buf_ != ''):
            if self.processor:
                #this.processor.process(buf_)
                buf_ += self.response.readline()
                [self.processor.process(line) for line in buf_.split('\n')]
            else:
                self.outfile_handle.write(buf_)
            buf_ = self.response.read(self.buffer_size)

    def fill_buffer(self):
        #extend parent version to catch HTTP specific error
        try:
            return SnaptronIterator.fill_buffer(self)
        except httplib.IncompleteRead, ir:
            sys.stderr.write(ir.partial)
            raise ir

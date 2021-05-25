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
import clsnapconf

class SnaptronIterator():

    def __init__(self,query_param_string,instance,endpoint,processor=None):
        self.buffer_size = clsnapconf.BUFFER_SIZE_BYTES
        self.query_param_string = query_param_string
        self.instance = instance
        self.endpoint = endpoint
        self.processor = processor

        self.idx = -1
        self.total_count = 0
        self.lines = []

        self.next = self.__next__

    def __iter__(self):
        return self
    
    def __next__(self):
        self.idx+=1
        lines_read = 0
        if self.idx >= len(self.lines):
            lines_read = self.fill_buffer()
            self.total_count += lines_read
            self.idx=0;
        if self.idx != 0 or lines_read > 0:
            if self.processor:
                return this.processor.process(self.lines[self.idx])
            return self.lines[self.idx]
        raise StopIteration

    def fill_buffer(self):
        buf_ = self.response.read(self.buffer_size)
        if buf_ is None or buf_ == '':
            return 0
        bufs = [buf_]
        last_char = buf_[-1]
        #top up to next newline
        while(last_char != None and last_char != '\n'):
            last_char = self.response.read(1)
            bufs.append(last_char)
        buf_ = "".join(bufs)
        #get rid of last newline so we dont get an emtpy string as a last array element
        if buf_[-1] == '\n':
            buf_ = buf_[:-1]
        self.lines = buf_.split("\n")
        return len(self.lines)

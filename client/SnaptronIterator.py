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
import clsnapconf

class SnaptronIterator():

    def __init__(self,query_param_strings,instances,endpoints):
        self.buffer_size = clsnapconf.BUFFER_SIZE_BYTES
        #multiple queries can either be 1) multiple junction ranges (e.g. both sides of an exon) or 2) multiple types: junction range + base range
        self.query_param_strings = query_param_strings
        self.instances = instances
        self.endpoints = endpoints

        self.total_count = 0
        self.query_idx = 0
        self.__init_for_query__()

        self.next = self.__next__
   
    def __init_for_query__(self):
        self.idx = -1
        self.lines = []

    def construct_query_string(self):
        pass
    
    def execute_query_string(self):
        pass

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
            return self.lines[self.idx]
        self.query_idx += 1
        if self.query_idx < len(self.query_param_strings):
            self.__init_for_query__()
            self.construct_query_string()
            self.execute_query_string()
            return self.__next__()
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
    

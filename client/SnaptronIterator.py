#!/usr/bin/env python

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
        if buf_ is None or buf_ == b'':
            return 0
        bufs = [buf_]
        last_char = buf_[-1]
        #top up to next newline
        while(last_char != None and last_char != b'\n'):
            last_char = self.response.read(1)
            bufs.append(last_char)
        buf_ = b''.join(bufs)
        #get rid of last newline so we dont get an emtpy string as a last array element
        if buf_[-1] == b'\n':
            buf_ = buf_[:-1]
        self.lines = buf_.split(b'\n')
        return len(self.lines)

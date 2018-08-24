#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 22:45:09 2018

@author: gustavo
"""

"""
Utility function to do a specific task.
"""

import re

def parse_GRDECL_DATA(fn):

    f = open(fn,'r')
    g = open('temp','w')    
    
    for line in f:
    
        # remove undesidered trailing hyphens found in COORDS lines in UNISIM file
        s = re.search('--\s\d+,\s\d+',line)
        if s is not None:
            line = line[0:s.start()] + '\n'
        g.write(line)
    
    
    f.close()
    g.close()
    
    
fn = '/Users/gustavo/projects/wellps/benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA'
parse_GRDECL_DATA(fn)
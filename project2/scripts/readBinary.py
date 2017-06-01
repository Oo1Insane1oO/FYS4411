#!/usr/bin/env python
# -*- coding: ascii -*-

import numpy as np
import sys

filename = sys.argv[1];
data = np.fromfile(filename, dtype=float, count=-1);

newData = np.zeros((len(data),6))

for i in range(len(data)/6):
    newData[i] = data[i:i+6]

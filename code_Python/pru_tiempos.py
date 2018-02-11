# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 09:59:13 2017

@author: algarra
"""

import numpy as np
import datetime
ciclos =10000000

start_time = datetime.datetime.now().time().strftime('%H:%M:%S')
for i in range(ciclos):
    a = np.random.uniform(0,0.5)
print("uniform")
end_time = datetime.datetime.now().time().strftime('%H:%M:%S')
total_time=(datetime.datetime.strptime(end_time,'%H:%M:%S') - datetime.datetime.strptime(start_time,'%H:%M:%S'))
print(total_time)
for i in range(ciclos):
    a = np.random.binomial(0,0.5)
print("binomial")
end_time = datetime.datetime.now().time().strftime('%H:%M:%S')
total_time=(datetime.datetime.strptime(end_time,'%H:%M:%S') - datetime.datetime.strptime(start_time,'%H:%M:%S'))
print(total_time)

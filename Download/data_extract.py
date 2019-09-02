# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from DynaData import *
import matplotlib.pyplot as plt
if __name__ == '__main__':
    PATH = r'C:\Users\chaos\Desktop\test'
    FILE = '\\nodout'
    NodoutFile = nodout_SMP(PATH+FILE,0,0,0,0)
    a = NodoutFile.simply
    FILE = '\\bndout'
    BndoutFile = bndout(PATH+FILE,0,0,0,0)
    b = BndoutFile.run
    bb = list(range(len(b['yforce'].tolist())))
    plt.plot(bb[10:],b['ytotal'][10:])
    plt.show()
    

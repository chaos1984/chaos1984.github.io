# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import Autolivlib as alv
if __name__ == '__main__':
    PATH = r'C:\Users\yujin.wang\Desktop\111'
    FILE = '\\nodout'
    NodoutFile = alv.nodout_SMP(PATH+FILE,0,0,0,0)
    NodoutFile.simply.to_csv("data.csv")
    FILE = '\\bndout'
    BndoutFile = alv.bndout(PATH+FILE,0,0,0,0)
    BndoutFile.run.to_csv("data1.csv")
    
    

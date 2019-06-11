# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

def fun(eps,K,n):
	return K*pow(eps,n)

data = pd.read_csv("EPSSTRAIN_TRUESTRESS.csv")
x = data.strain
y = data.stress - data.stress[0]
#plt.figure(1)
#plt.xlabel("strain")
#plt.ylabel("stress")
#plt.tittl("Efecctive Plastic Strain Vs. True Stress")
p0 = 10000,0.4
popt,pcov = sp.optimize.curve_fit(fun, x, y, p0)
K,n = popt
x_fit = np.linspace(0,0.22,100)
plt.plot(x,y,'r')
plt.plot(x_fit,fun(x_fit,K,n),'b')
plt.show()
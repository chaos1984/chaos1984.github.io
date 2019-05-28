# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Sample information
l = 10. #Gauge length
h = 3   #Thickness
tb = 10.#Width
A = h*tb
curveid = 2350

test_data = pd.read_csv("data.csv")

COL = len(test_data.columns) 

fout_key = open("Curve.key",'w')
fout_fd = open("fd_test.txt","w")
listE = [] 
listFailureStrain = []
plt.figure()
fout_key.write("*KEYWORD\n*PARAMETER\n")
fout_key.write("*DEFINE_TABLE\n%d\nINPUTYOURRATE\n" %(curveid))
for i in range(COL/2):
	
	ES = test_data[test_data.columns[2*i+1]].dropna()/1000
	EE = test_data[test_data.columns[2*i]].dropna()/10
	F = ES*A
	D = EE*l
	TE = np.log(1+EE)
	
	TS = ES*(1+EE)
	
	listE.append(TS[1]/TE[1])
	
	EPS = TE- TS/listE[i]
	listFailureStrain.append(EPS.tolist()[-1])
	
	plt.subplot(211)
	plt.plot(EPS,TS)
	plt.legend([100,10,1])
	plt.subplot(212)
	plt.plot(D,F)
	plt.legend([100,10,1])
	curveid+=1
	fout_key.write("*DEFINE_CURVE_TITLE\n%s\n%d,0,1.0,1.0,0.0,0.0,0\n" %(test_data.columns[2*i+1],curveid))
	for i in range(len(EPS)-1):
		fout_key.write("%f,%f\n" %(EPS[i+1],TS[i+1]))
		fout_fd.write("%f,%f\n" %(D[i+1],F[i+1]))
print np.mean(listE)
fout_key.write("*DEFINE_CURVE_TITLE\nStrainRate VS. FailureStrain\n%d,0,1.0,1.0,0.0,0.0,0\n" %(curveid+1))
for i in listFailureStrain:
	fout_key.write("Rate,%f\n" %(i))
fout_key.write("*DEFINE_CURVE_TITLE\nStressTriaxiality VS. FailureStrain\n%d,0,1.0,1.0,0.0,0.0,0\n" %(curveid+2))
for i in listFailureStrain:
	fout_key.write("-1,1\n")
fout_key.write("*END")
plt.show()
fout_key.close()
fout_fd.close()
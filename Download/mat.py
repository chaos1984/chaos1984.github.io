# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
import scipy.signal as signal

#def fun(eps,a0,a1,a2,a3):
#	return a0*(pow(eps,a1)+eps*a2+pow(eps,a3))
#def fun(eps,a0,a1,a2,a3,a4):
#	return a0*np.arctan(a1*eps) + pow(eps,a2) 
def fun(eps,a0,a1,a2,a3,a4):
	return a0 + a1*pow(eps,1) + a2*pow(eps,2)+a3*pow(eps,3)+a4*pow(eps,4)

def MAT89LCSS(fdfile,A,l0):
    data = pd.read_csv(fdfile)
    x = data.Disp
    y = data.Force
    EngStrain = x/l0
    EngStress = y/A
    EngStressFilt =  pd.Series(signal.medfilt(EngStress,201))
    PeakStress = max(EngStressFilt)
    PeakStrain = EngStrain[EngStressFilt.idxmax()]
    print "PeakStress:%f\tPeakStrain:%f\t Elastic modulus:%f" %(PeakStress,PeakStrain,PeakStress/PeakStrain)
    plt.figure(1)
    plt.grid("on")
    plt.title("Engineering Strain VS. Engineering Stress Curve")
    plt.xlabel("Engineering Strain[-]")
    plt.ylabel("Engineering Stress[GPa]")
#    plt.plot(PeakStrain,PeakStress,'p')
    plt.plot(EngStrain,EngStress)
#    plt.plot(EngStrain,EngStressFilt,'--')
    plt.legend(["Strain Rate 0.001/ms","Strain Rate 0.01/ms","Strain Rate 0.1/ms","Strain Rate 1/ms"])
    TrueStrain =  np.log(1+EngStrain)
    TrueStress = EngStressFilt*(1+EngStrain)
    plt.figure(2)
#    plt.plot(EngStrain,EngStress,'b')
    plt.grid("on")
    plt.title("True Strain VS. True Stress Curve")
    plt.xlabel("True Strain[-]")
    plt.ylabel("True Stress[GPa]")
    plt.plot(TrueStrain,TrueStress)
    plt.legend(["Strain Rate 0.001/ms","Strain Rate 0.01/ms","Strain Rate 0.1/ms","Strain Rate 1/ms"])
    return TrueStrain,TrueStress,[PeakStrain,PeakStress]

def MAT24LCSS(x,y,p,length=2.0,yieldpoint=0.02,Ymould_user=0,Peak=[0,0],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0,alignstrain=1.45,pointnum=50):
    Ymould = Ymould_user
    print " User defined Ymould:",Ymould
    y_eps = []
    x_eps = []
    for index,value in enumerate(x):
        if value >= yieldpoint:
            x_eps.append(value - y[index]/Ymould)
            y_eps.append(y[index])
    plt.figure(3)
    displayratio = int(len(x_eps)*ratio)
    plt.plot(x_eps[:displayratio],y_eps[:displayratio],'--')
    x_eps_mod = x_eps;y_eps_mod = y_eps
    delta = int(len(x_eps)/pointnum)
    x_eps_mod = [value for index,value in enumerate(x_eps) if index%delta==0];y_eps_mod = [value for index,value in enumerate(y_eps) if index%delta==0]
    for index,value in enumerate(x_eps_mod):
        if value > strainturn:
            y_eps_mod[index]= y_eps_mod[index-1]+scaleturn*(x_eps_mod[index]-x_eps_mod[index-1])
        else:
            continue
    if x_eps_mod[-1] < alignstrain:
        x_eps_mod.append(alignstrain)
        y_eps_mod.append(y_eps_mod[-1]+scaleturn*(x_eps_mod[-1]-x_eps_mod[-2]))
    for index,value in enumerate(x_eps_mod):
        y_eps_mod[index] = y_eps_mod[index]*curvescale
    for index in range(len(y_eps_mod)):
        try:
            value = y_eps_mod[index]
            if value > y_eps_mod[index+1]:
               y_eps_mod.pop(index)
               x_eps_mod.pop(index)
               index -= 1
        except:
            pass
    if extend != 0:
        x_eps_mod.append(x_eps_mod[-1]+extend)
        y_eps_mod.append(y_eps_mod[-1]+extend*(scaleextend))
    plt.title("Modified Effective Plastic Strain VS. Stress Curve")
    plt.xlabel("Effective Plastic Strain[-]")
    plt.ylabel("Stress[GPa]")
    plt.plot(x_eps_mod[:displayratio],y_eps_mod[:displayratio])
    plt.legend(["Strain Rate 0.001/ms","Modified Strain Rate 0.001/ms","Strain Rate 0.01/ms","Modified  Strain Rate 0.01/ms","Strain Rate 0.1/ms","Modified Strain Rate 0.1/ms","Strain Rate 1/ms","Modified Strain Rate 1/ms",])
    return x_eps_mod,y_eps_mod


if __name__ == '__main__':
    Ymould_user = 1.5
    yieldpoint= 0.1
    pointnum = 80
    strainturnlist=[1.2,1.2,1.2,0.75];scaleturnlist=[0.09,0.09,0.085,0.060];curvescalelist = [1,1.,1.03,1.03]
    A = 30; l0 = 10
    FDfile = ["test10.txt","test100.txt","test1000.txt","test10000.txt"]
    strainrange = 1.3;extend = 0.5;scaleextend=0.08;alignstrain =1.45
    ratio = 1
    fout = open("Curve.key",'w')
    strain_rate = [-6.907,-4.605,-2.302,0]
    curvenum = 2350
    fout.write('*KEYWORD\n')
    fout.write('$ Created: ' + time.strftime("%d.%m.%Y %H:%M:%S") + '\n')
    fout.write('$ Parameters:\n$ A:%f\n$ l0:%f\n$ Fitting range:%f\n$ extend strain:%f\n$ scaleextend:%f\n$ alignstrain:%f\n' %(A,l0,strainrange,extend,scaleextend,alignstrain))
    fout.write("$ Control List:\n$ strain turning:%s\n$ strain scaleturn:%s\n$ curvescale:%s\n" %(str(strainturnlist),str(scaleturnlist),str(curvescalelist)))
    fout.write("*DEFINE_TABLE\n%d\n" %(curvenum))
    for i in strain_rate:
        fout.write("%f\n" %(i))    
    for index,File in enumerate(FDfile):
        curvenum += 1
        strainturn = strainturnlist[index];scaleturn = scaleturnlist[index];curvescale = curvescalelist[index]
        x,y,peak = MAT89LCSS(File,A,l0) 
        x_fit,y_fit= MAT24LCSS(x,y,p0,strainrange,yieldpoint=yieldpoint,Ymould_user=Ymould_user,Peak=peak,strainturn=strainturn,scaleturn=scaleturn,curvescale=curvescale,ratio=ratio,extend=extend,scaleextend=scaleextend,alignstrain =alignstrain,pointnum=pointnum)    
        fout.write('*DEFINE_CURVE_TITLE\nRate %.5f\t%s\n' %(pow(np.e,strain_rate[index]),File))
        fout.write('$     LCID      SIDR       SFA       SFO      OFFA      OFFO    DATTYP\n')
        fout.write('      %d         0    1.0000&scale        0.0000    0.0000\n' %(curvenum))
        for i in range(len(x_fit)):
            fout.write("%f,%f\n" %(x_fit[i],y_fit[i]))
    fout.write("*END\n")
    fout.close()
#    plt.plot(x_fit,y_fit-a3,'y')# Fitting cuver for dyna MAT24    
    plt.grid('on')
    plt.show()
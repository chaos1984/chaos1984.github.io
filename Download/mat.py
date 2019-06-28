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
    PeakStress = max(EngStress)
    PeakStrain = EngStrain[EngStress.idxmax()]
    print "PeakStress:%f\tPeakStrain:%f" %(PeakStress,PeakStrain)
    plt.figure(5)
    plt.grid("on")
    plt.title("Engineering Strain VS. Engineering Stress Curve")
    plt.xlabel("Engineering Strain[-]")
    plt.ylabel("Engineering Stress[GPa]")
    plt.plot(EngStrain,EngStress)
    TrueStrain = np.log(1+EngStrain)
    TrueStress = EngStress*(1+EngStrain)
    plt.figure(1)
#    plt.plot(EngStrain,EngStress,'b')
    plt.plot(TrueStrain,TrueStress,"--")
    return TrueStrain,TrueStress,PeakStress,PeakStrain

def MAT24LCSS(x,y,PeakStress,PeakStrain,p,length=2.0,yieldfpoint=0.02,Ymould_user=0,scale=[1,1,1,1,1],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0):
    popt,pcov = sp.optimize.curve_fit(fun, x, y, p0)
    a0,a1,a2,a3,a4 = popt
    SigmaY = fun(yieldfpoint,a0,a1,a2,a3,a4)
    if Ymould_user != 0:
        Ymould = Ymould_user
        print " User defined Ymould:",Ymould
    else:
        Ymould = SigmaY/yieldfpoint
        print "Yield stress:",SigmaY,"Ymould:",Ymould
    x_fit = np.linspace(0,length,500)
    y_fit = fun(x_fit,a0,a1,a2,a3,a4)
    plt.plot(x_fit,y_fit,'r')#fitting cuve for test
    y_fit_eps = []
    x_fit_eps = []
    for i in x_fit:
        if i >= 0.2:
            y_eps = fun(i,a0,a1,a2,a3,a4)
            x_eps = i - y_eps/Ymould
            x_fit_eps.append(x_eps)
            y_fit_eps.append(y_eps)
    popt,pcov = sp.optimize.curve_fit(fun,x_fit_eps,y_fit_eps, p0)
    a0,a1,a2,a3,a4 = popt
    print ("a0:%8.5f,a1:%5.3f,a2:%5.3f,a3:%5.3f,a4:%5.3f" %(a0,a1,a2,a3,a4))
    x_fit = np.linspace(0,length+extend,50)
    y_fit = fun(x_fit,a0,a1,a2,a3,a4)
    y_fit_modified = fun(x_fit,a0*scale[0],a1*scale[1],a2*scale[2],a3*scale[3],a4*scale[4])
    for index,value in enumerate(x_fit):
        if value < np.log(1+PeakStrain):
            y_fitmodified
        if value > strainturn:
            y_fit_modified [index]= y_fit_modified[index-1]+scaleturn*(x_fit[index]-x_fit[index-1])
        else:
            continue
        if value > length:
            y_fit_modified [index]= y_fit_modified[index-1]+scaleextend*(x_fit[index]-x_fit[index-1])
        else:
            continue
    for index,value in enumerate(x_fit):
        y_fit_modified [index] = y_fit_modified [index]*curvescale
    plt.figure(2)
    plt.plot(x_fit_eps,y_fit_eps,'b')#MAT24
    plt.figure(3)
    displayratio = int(len(x_fit)*ratio)
    plt.title("Modified Effective Plastic Strain VS. Stress Curve")
    plt.xlabel("Effective Plastic Strain[-]")
    plt.ylabel("Stress[GPa]")
    plt.plot(x_fit[:displayratio],y_fit[:displayratio],"--")
    plt.plot(x_fit[:displayratio],y_fit_modified[:displayratio])# Fitting cuver for dyna MAT24
    return x_fit,y_fit_modified,popt


if __name__ == '__main__':
    Ymould_user = 1.5
    yieldfpoint= 0.02
    scale = [1,1,1,1,1]
    strainturnlist=[1.2,1.2,1.2,0.75];scaleturnlist=[0.09,0.09,0.085,0.060];curvescalelist = [1,1.,1.012,1.]
    A = 30; l0 = 10
    FDfile = ["test10.txt","test100.txt","test1000.txt","test10000.txt"]
    strainrange = 1.45;extend = 1.0;scaleextend=0.08
    ratio = 0.6
    fout = open("Curve.key",'w')
    strain_rate = [-6.907,-4.605,-2.302,0]
    curvenum = 2350
    p0 = 1,1,1,1,1
    fout.write('*KEYWORD\n')
    fout.write('$ Created: ' + time.strftime("%d.%m.%Y %H:%M:%S") + '\n')
    fout.write('$ Parameters:\n$ A:%f\n$ l0:%f\n$ yieldfpoint:%f\n$ Fitting range:%f\n$ extend strain:%f\n$ scaleextend:%f\n' %(A,l0,yieldfpoint,strainrange,extend,scaleextend))
    fout.write("$ Control List:\n$ strain turning:%s\n$ strain scaleturn:%s\n$ fitting function coefficents scale:%s\n$ curvescale:%s\n" %(str(strainturnlist),str(scaleturnlist),str(scale),str(curvescalelist)))
    fout.write("*DEFINE_TABLE\n%d\n" %(curvenum))
    for i in strain_rate:
        fout.write("%f\n" %(i))    
    for index,File in enumerate(FDfile):
        curvenum += 1
        strainturn = strainturnlist[index];scaleturn = scaleturnlist[index];curvescale = curvescalelist[index]
        x,y,PeakStress,PeakStrain = MAT89LCSS(File,A,l0) 
        x_fit,y_fit,a = MAT24LCSS(x,y,PeakStress,PeakStrain,p0,strainrange,yieldfpoint=yieldfpoint,Ymould_user=Ymould_user,scale=scale,strainturn=strainturn,scaleturn=scaleturn,curvescale=curvescale,ratio=ratio,extend=extend,scaleextend=scaleextend)    
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
# -*- coding: utf-8 -*-
"""
Author:Yujin.Wang
Date:2019.04.07
The lib. is used to translate the F-D curve to MAT24 effective plastic strain - true stress.
Usage:
MAT89LCSS(fdfile,A,l0) 
return TrueStrain,TrueStress,[PeakStrain,PeakStress]
MAT24LCSS(truestrain,truestress,yieldpoint=0.02,Ymould_user=0,Peak=[0,0],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0,alignstrain=1.45,pointnum=50):
return x_eps_mod,y_eps_mod
CurveKey(A,l0,FDfile,strain_rate,curvenum,ratio,Ymould_user,yieldpoint,pointnum,strainturnlist,scaleturnlist,curvescalelist,alignstrain,extend,scaleextend):
No return
################################################################################
Parameters:
    A - float, Specemen Crosssection area
	l0 - float, Gage lengrh
    FDfile - Filename, Force-Disp curve
    strain_rate - list, Strain rate
    curvenum - int, Curve number in .key file
    ratio - float, Display ratio
#################################User Define###################################
    Ymould_user - float, Young modulus
    yieldpoint - float, Yield point
    pointnum - int, the number of output points
    strainturnlist - list, Turn strain 
    scaleturnlist - list,  Scale of the turn strain 
    curvescalelist - list, Curve scale equalent to SFO
    alignstrain - float, align the strain to user defined value(>max(strainturnlist))
    extend - float, Extend strain to a user defined value
    scaleextend - float, Scale of the extend strain 
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate,optimize
import scipy.signal as signal
from sympy import symbols,nsolve

def f_eps(eps,c1,c2,a):
    return (np.power(np.e,-c1*eps) + np.power(eps,c2))*(1-np.power(np.e,-a*eps))

def h_epsrate_T(epsrate,m,T,alpha):
    return np.power(epsrate,m)*np.power(np.e,alpha/T)

def sigma_eps_epsrate_T(eps,epsrate,T,K,c1,c2,c3,c4,a,m,alpha):
    f = f_eps(eps,c1,c2,a)
    h = h_epsrate_T(epsrate,m,T,alpha)
    value = K * ( f + (eps*np.power(np.e,(1-eps/(c3*h)))/(c3*h)-f) * np.power(np.e,(np.log(h)-c4)*eps))*h
#    print (eps*np.power(np.e,(1-eps/(c3*h)))/(c3*h)-f)
    return value


def DGSZ_c1c2(eps,sigma):
    c1,c2 = symbols('c1,c2')
    equations = [(np.power(np.e,-c1*eps[0]) + np.power(eps[0],c2))     \
                 /(np.power(np.e,-c1*eps[1]) + np.power(eps[1],c2)) - sigma[0]/sigma[1], \
                 (np.power(np.e,-c1*eps[0]) + np.power(eps[0],c2))     \
                 /(np.power(np.e,-c1*eps[1]) + np.power(eps[2],c2)) - sigma[0]/sigma[2]]
    result = nsolve(equations,[c1,c2],[0,0])
    return result
    
def DGSZ_m(epsrate,sigma):
    return np.log(sigma[0]/sigma[1]) / np.log(epsrate[0]/epsrate[1])

def DGSZ_alpha(T,sigma):
    return np.log(sigma[0]/sigma[1])/(1./T[0]-1/T[1])

def DGSZ_K(sigma,eps,epsrate,m,T,alpha,c1,c2):
    h = h_epsrate_T(epsrate,m,T,alpha)
    return sigma/(h*np.power(np.e,-c1*eps)+np.power(eps,-c2))

def DGSZ_c3(eps,epsrate,m,T,alpha):
    h = h_epsrate_T(epsrate,m,T,alpha)
    return eps/h
    
def DGSZ_c4(epsrate,m,alpha,T):
    return 200. + np.log(np.power(epsrate,m)*np.power(np.e,alpha/T))

def DSGZ_a(eps):
    return -np.log(0.03)/eps

def CowperSymondsFunc(x,c,p):
    rate,ES0 = x
    return ES0*(1+np.power((rate/c),1./p))


def MAT89LCSS(fdfile,A,l0):
    data = pd.read_csv(fdfile)
    x = data.Disp
    y = data.Force
    EngStrain = x/l0
    EngStress = y/A
    EngStressFilt =  pd.Series(signal.medfilt(EngStress,201))
    PeakStress = max(EngStressFilt[:800])
    PeakStrain = EngStrain[EngStressFilt[:800].idxmax()]
    print "PeakStress:%f\tPeakStrain:%f\t Elastic modulus:%f" %(PeakStress,PeakStrain,PeakStress/PeakStrain)
    plt.figure(1)
    plt.grid("on")
    plt.scatter(PeakStrain,PeakStress, marker='o')
    plt.title("Engineering Strain VS. Engineering Stress Curve")
    plt.xlabel("Engineering Strain[-]")
    plt.ylabel("Engineering Stress[GPa]")
#    plt.xticks([])
#    plt.yticks([])
#    plt.plot(PeakStrain,PeakStress,'p')
    plt.plot(EngStrain,EngStressFilt)
#    plt.plot(EngStrain,EngStressFilt,'--')
    plt.legend(["0.001mm/ms","0.01mm/ms","0.1mm/ms","1mm/ms"])
    TrueStrain =  np.log(1+EngStrain)
    TrueStress = EngStressFilt*(1+EngStrain)
    plt.figure(2)
#    plt.plot(EngStrain,EngStress,'b')
    plt.grid("on")
    plt.title("True Strain VS. True Stress Curve")
#    plt.xticks([])
#    plt.yticks([])
    plt.xlabel("True Strain[-]")
    plt.ylabel("True Stress[GPa]")
    plt.plot(TrueStrain,TrueStress)
    plt.legend(["0.001mm/ms","0.01mm/ms","0.1mm/ms","1mm/ms"])
    return TrueStrain,TrueStress,[PeakStrain,PeakStress]

def MAT24LCSS(x,y,yieldpoint=0.02,Ymould_user=0,Peak=[0,0],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0,alignstrain=1.45,pointnum=50,smoothflag=0):
    print " User defined Ymould:",Ymould_user
    y_eps = []
    x_eps = []
    for index,value in enumerate(x):
        if value >= yieldpoint and (value > x[index-1] or y[index]>y[index-1]):
            x_eps.append(value - y[index]/Ymould_user)
            y_eps.append(y[index])
    plt.figure(4)
    displayratio = int(len(x_eps)*ratio)
#    plt.plot(x_eps[:displayratio],y_eps[:displayratio],'--')
    x_eps_mod =[0] 
    y_eps_mod = [y_eps[0]/3.]
    delta = int(len(x_eps)/pointnum)
    x_eps_mod.extend( [value for index,value in enumerate(x_eps) if index%delta==0]);y_eps_mod.extend([value for index,value in enumerate(y_eps) if index%delta==0])
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
    if extend != 0:
        x_eps_mod.append(x_eps_mod[-1]+extend)
        y_eps_mod.append(y_eps_mod[-1]+extend*(scaleextend))
        plt.title("Modified Effective Plastic Strain VS. Stress Curve")
    plt.xlabel("Effective Plastic Strain[-]")
    plt.ylabel("Stress[GPa]")
    plt.grid('on')
    if smoothflag == 1:
        print "The curve has been smoothed!"
        f = interpolate.interp1d(x_eps_mod,y_eps_mod,kind='cubic')
        x_eps_smooth = np.linspace(x_eps_mod[0],x_eps_mod[-1],pointnum)
        y_eps_smooth = f(x_eps_smooth)
        plt.plot(x_eps_smooth[:displayratio],y_eps_smooth[:displayratio])
#        plt.legend(["Strain Rate 0.001/ms","Modified Strain Rate 0.001/ms","Strain Rate 0.01/ms","Modified  Strain Rate 0.01/ms","Strain Rate 0.1/ms","Modified Strain Rate 0.1/ms","Strain Rate 1/ms","Modified Strain Rate 1/ms",])
        return x_eps_smooth,y_eps_smooth
    elif smoothflag == 0:
        plt.plot(x_eps_mod[:displayratio],y_eps_mod[:displayratio])
#        plt.legend(["Strain Rate 0.001/ms","Modified Strain Rate 0.001/ms","Strain Rate 0.01/ms","Modified  Strain Rate 0.01/ms","Strain Rate 0.1/ms","Modified Strain Rate 0.1/ms","Strain Rate 1/ms","Modified Strain Rate 1/ms",])
        return x_eps_mod,y_eps_mod

def CurveKey(A,l0,FDfile,strain_rate,curvenum,ratio,Ymould_user,yieldpoint,pointnum,strainturnlist,scaleturnlist,curvescalelist,alignstrain,extend,scaleextend,smooth=0):
    fout = open("Curve.key",'w')
    fout.write('*KEYWORD\n')
    fout.write('$ Created: ' + time.strftime("%d.%m.%Y %H:%M:%S") + '\n')
    fout.write('$ Parameters:\n$ A:%f\n$ l0:%f\n$ Young modulus:%f\n$ Yield point:%f\n$ Alignstrain:%f\n$ Extend strain:%f\n$ Scaleextend:%f\n' %(A,l0,Ymould_user,yieldpoint,alignstrain,extend,scaleextend))
    fout.write("$ Control List:\n$ Strain turning:%s\n$ Strain scaleturn:%s\n$ Curvescale:%s\n" %(str(strainturnlist),str(scaleturnlist),str(curvescalelist)))
    fout.write("*DEFINE_TABLE\n%d\n" %(curvenum))
    for i in strain_rate:
        fout.write("%f\n" %(i))    
    for index,File in enumerate(FDfile):
        curvenum += 1
        strainturn = strainturnlist[index];scaleturn = scaleturnlist[index];curvescale = curvescalelist[index]
        x,y,peak = MAT89LCSS(File,A,l0) 
        x_fit,y_fit= MAT24LCSS(x,y,yieldpoint=yieldpoint,Ymould_user=Ymould_user,Peak=peak,strainturn=strainturn,scaleturn=scaleturn,curvescale=curvescale,ratio=ratio,extend=extend,scaleextend=scaleextend,alignstrain =alignstrain,pointnum=pointnum,smoothflag=smooth)    
        fout.write('*DEFINE_CURVE_TITLE\nRate %.5f\t%s\n' %(pow(np.e,strain_rate[index]),File))
        fout.write('$     LCID      SIDR       SFA       SFO      OFFA      OFFO    DATTYP\n')
        fout.write('      %d         0    1.0000&scale        0.0000    0.0000\n' %(curvenum))
        for i in range(len(x_fit)):
            fout.write("%f,%f\n" %(x_fit[i],y_fit[i]))
    fout.write("*END\n")
    fout.close()
#    return x_fit,y_fit
#    plt.show()

#def CowperSymondsCurve(File,A,l0,Ymould_user,pointnum,ratelist,x_p,y_p):
#    y_eps = []
#    x_eps = []
#    for index,value in enumerate(x):
#        if value >= yieldpoint and (value > x[index-1] or y[index]>y[index-1]):
#            x_eps.append(value - y[index]/Ymould_user)
#            y_eps.append(y[index])
#    plt.figure(4)
#    plt.plot(x_eps,y_eps,'--')
#    popt,pcov = optimize.curve_fit(CowperSymondsFunc,x_p,y_p,[0.001,0.001])
#    c = popt[0];p = popt[1]
#    print 'c:',c
#    print 'p:',p
##    print "Cowper-Symonds Constants:\nc:%f\tp:%f\n" %(c,p)
#    x_eps_fit = np.linspace(0,x_eps[-1],pointnum)
#    for rate in ratelist:
#        plt.figure(4)
#        y_eps_fit = CowperSymondsFunc(x_eps_fit,c,p,rate=0.001)
#        plt.plot(x_eps_fit,y_eps_fit)
#    return x_eps_fit,y_eps_fit
    

if __name__ == '__main__':
    A = 30; l0 = 10 #样件截面积，标距
    FDfile = ["test10.txt","test100.txt","test1000.txt","test10000.txt"]#力位移曲线
    strain_rate = [-6.907,-4.605,-2.302,0] #力位移曲线对应的应变率
    curvenum = 2350 #key文件中曲线的编号起始编号
    ratio = 1 #显示比例
#################################User Define###################################
    Ymould_user = 1.5 #定义弹性模量
    yieldpoint= 0.105 #屈服点（避开交叉区域）
    pointnum = 80     #输出点的个数
    strainturnlist=[1.3,1.3,1.35,0.9]#指定每根曲线发生转折的应变
    scaleturnlist=[0.08,0.08,0.06,0.0450];#指定每根曲线发生转折的比例
    curvescalelist = [1,1.05,1.05,1.1]#指定每个区县整体偏移的比例
#    strainturnlist=[1.2,1.2,1.2,0.75];scaleturnlist=[0.06,0.06,0.055,0.050];curvescalelist = [1,1.,1.05,1.05]
    alignstrain = 1.45#应力应变曲线对齐的位置（应变值）
    extend = 0.5     #对齐后曲线延长距离（应变值）
    scaleextend=0.08 #延长段的增量
#################################User Define###################################
    CurveKey(A,l0,FDfile,strain_rate,curvenum,ratio,Ymould_user,yieldpoint,pointnum,strainturnlist,scaleturnlist,curvescalelist,alignstrain,extend,scaleextend)
    ratelist = [np.power(np.e,i) for i in strain_rate]
    y_p = []
    x_eps=[];y_eps=[]
    for index,File in enumerate(FDfile):
        x,y,peak = MAT89LCSS(File,A,l0)
        y_p.append(peak[1])
        if index == 0:
            x0 = x
            y0 = y
    for index,value in enumerate(x0):
        if value >= yieldpoint and (value > x0[index-1] or y0[index]>y[index-1]):
            x_eps.append(value - y0[index]/Ymould_user)
            y_eps.append(y0[index])
    popt,pcov = optimize.curve_fit(CowperSymondsFunc,[ratelist,y_eps[0]],y_p,[0.1,0.1])
    c,p = popt
    print 'c:',c
    print 'p:',p
    plt.figure(5)
    plt.scatter(ratelist,y_p,color='r',marker='*')
    y_p_fit=CowperSymondsFunc([ratelist,y_p[0]],c,p)
    plt.scatter(ratelist,y_p_fit)
#    x_eps_fit = np.linspace(0,x.tolist()[-1],pointnum)
    plt.figure( 3)
    fout = open("Curve.key",'w')
    fout.write('*KEYWORD\n')
    fout.write("*DEFINE_TABLE\n%d\n" %(curvenum))
    for i in strain_rate:
        fout.write("%f\n" %(i))    
    for rate in ratelist:
        curvenum += 1
        y_eps_fit = CowperSymondsFunc([rate,np.array(y_eps)],c,p)
        fout.write('*DEFINE_CURVE_TITLE\nRate %.5f\n' %(rate))
        fout.write('$     LCID      SIDR       SFA       SFO      OFFA      OFFO    DATTYP\n')
        fout.write('      %d         0    1.0000&   1.0000    0.0000    0.0000\n' %(curvenum))
        for i in range(len(x_eps)):
            fout.write("%f,%f\n" %(x_eps[i],y_eps_fit[i]))
        plt.plot(x_eps,y_eps_fit,dashes=[2, 2, 10, 2])
    fout.write("*END\n")
    
    
    ###########################################################################
    eps = [0.2,0.4,0.6];sigma=[47.7,55.5,62.2];T = 238.
    c1,c2 = DGSZ_c1c2(eps,sigma)
    epsrate = [1000.,100.];sigma =[47.7,36.7]
    m = DGSZ_m(epsrate,sigma)
    T = [238.,238.1];sigma=[47.7,47.6]
    alpha = DGSZ_alpha(T,sigma)
    T = T[0];epsrate = epsrate [0];sigma = 47.7;eps=0.2
    K = DGSZ_K(sigma,eps,epsrate,m,T,alpha,c1,c2)
    c3 = DGSZ_c3(eps,epsrate,m,T,alpha)
    c4 = DGSZ_c4(epsrate,m,alpha,T)
    a = DSGZ_a(0.2)
    ##########################
#    m_list =[0.05,0.075,0.1]
#    m = 0.075
#    c2= 1.
#    c1 = 0.5
#    c1_list = [0.2]
#    alpha_list = [800.,1000.,1200.]
#    K_list = [0.55]
#    c3_list = [0.001,0.1]
#    c4_list = [11.,15.,20.]
    ##########################
    plt.figure(4)
    epsrate_list = [10.,100.,1000.]
    T_list = [238.]
    for T in T_list:
        for epsrate in epsrate_list:
#            for K in K_list:
                eps = np.linspace(0,2,150)
                sigma = sigma_eps_epsrate_T(eps,epsrate,T,K,c1,c2,c3,c4,a,m,alpha)/1000
                y_eps = []
                x_eps = []
                for index,value in enumerate(eps):
                        x_eps.append(value - sigma[index]/1.5e3)
                plt.plot(x_eps,sigma,marker="*")
    plt.grid("on")
    print "c1\tc2\tm\ta\tK\tc3\tc4\talpha"
    print "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f" %(c1,c2,m,a,K,c3,c4,alpha)
    plt.show()
#    x,y = CowperSymondsCurve(FDfile[0],A,l0,Ymould_user,80,ratelist,x_p,y_p)
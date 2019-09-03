# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import Autolivlib as alv
import matplotlib.pyplot as plt
# -*- coding: utf-8 -*-
"""
Author:Yujin.Wang
Date:2019.04.07
The lib. is used to translate the F-D curve to MAT24 effective plastic strain - true stress.
Usage:
MAT89LCSS(fdfile,A,l0) 
return TrueStrain,TrueStress,[PeakStrain,PeakStress]
MAT24LCSS(truestrain,truestress,yieldpoint=0.02,Ymould_user=0,Peak=[0,0],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0,alignstrain=1.45,pointnum=50):
return x_strain_mod,y_strain_mod
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

def f_strain(strain,c1,c2,a):
    return (np.power(np.e,-c1*strain) + np.power(strain,c2))*(1-np.power(np.e,-a*strain))

def h_strainrate_T(strainrate,m,T,alpha):
    return np.power(strainrate,m)*np.power(np.e,alpha/T)

def stress_strain_strainrate_T(strain,strainrate,T,K,c1,c2,c3,c4,a,m,alpha):
    f = f_strain(strain,c1,c2,a)
    h = h_strainrate_T(strainrate,m,T,alpha)
    value = K * ( f + (strain*np.power(np.e,(1-strain/(c3*h)))/(c3*h)-f) * np.power(np.e,(np.log(h)-c4)*strain))*h
#    print (strain*np.power(np.e,(1-strain/(c3*h)))/(c3*h)-f)
    return value


def DGSZ_c1c2(strain,stress):
    c1,c2 = symbols('c1,c2')
    equations = [(np.power(np.e,-c1*strain[0]) + np.power(strain[0],c2))     \
                 /(np.power(np.e,-c1*strain[1]) + np.power(strain[1],c2)) - stress[0]/stress[1], \
                 (np.power(np.e,-c1*strain[0]) + np.power(strain[0],c2))     \
                 /(np.power(np.e,-c1*strain[1]) + np.power(strain[2],c2)) - stress[0]/stress[2]]
    result = nsolve(equations,[c1,c2],[0,0])
    return result
    
def DGSZ_m(strainrate,stress):
    return np.log(stress[0]/stress[1]) / np.log(strainrate[0]/strainrate[1])

def DGSZ_alpha(T,stress):
    return np.log(stress[0]/stress[1])/(1./T[0]-1/T[1])

def DGSZ_K(stress,strain,strainrate,m,T,alpha,c1,c2):
    h = h_strainrate_T(strainrate,m,T,alpha)
    return stress/(h*np.power(np.e,-c1*strain)+np.power(strain,-c2))

def DGSZ_c3(strain,strainrate,m,T,alpha):
    h = h_strainrate_T(strainrate,m,T,alpha)
    return strain/h
    
def DGSZ_c4(strainrate,m,alpha,T):
    return 200. + np.log(np.power(strainrate,m)*np.power(np.e,alpha/T))

def DSGZ_a(strain):
    return -np.log(0.03)/strain

def CowperSymondsFunc(x,c,p):
    rate,ES0 = x
    return ES0*(1+np.power((rate/c),1./p))

def ConditionRemove(data,condition):
    a =[]
    isflag = 1
    string ="if "+ condition +":\n"
    string +="                     for k in range(len(data)):\n   \
                        data[k].pop(i+1)\n   \
                  isflag = 1\n"
    while isflag == 1:
        a.extend(data)
        len0 = len(a[0])
        try:
            for i in range(len(a[0])):
                exec(string)
        except:
            if len(data[0]) == len0:
                break
    return data
    
    
def MAT89LCSS(fdfile,A,l0):
    '''fdfile -- Force-Displacement Curve
       A -- Area of the Cross section of specimen
       l -- Length of the specimen (gage length or grip length or equalent length)
       Output : TrueStrain,TrueStress,[PeakStrain,PeakStress]
     '''
    data = pd.read_csv(fdfile)
    x = data.Disp
    y = data.Force
    y=  pd.Series(signal.medfilt(y,201))
    print np.trapz(y,x)
    plt.figure(1)
    plt.grid("on")
    plt.title("Force VS. Displacement Curve")
    plt.xlabel("Displacement[mm]")
    plt.ylabel("Force[kN]")
    plt.plot(x,y)
    EngStrain = x/l0
    EngStress = y/A
    # EngStressFilt =  pd.Series(signal.medfilt(EngStress,201))
    PeakStress = max(EngStress[:800])
    PeakStrain = EngStrain[EngStress[:800].idxmax()]
    print "PeakStress:%f\tPeakStrain:%f\t Elastic modulus:%f" %(PeakStress,PeakStrain,PeakStress/PeakStrain)
    plt.figure(2)
    plt.grid("on")
    plt.scatter(PeakStrain,PeakStress, marker='o')
    plt.title("Engineering Strain VS. Engineering Stress Curve")
    plt.xlabel("Engineering Strain[-]")
    plt.ylabel("Engineering Stress[GPa]")
    plt.plot(EngStrain,EngStress)
    plt.legend(["0.01mm/ms","0.1mm/ms","1mm/ms","10mm/ms"])
    TrueStrain =  np.log(1+EngStrain)
    TrueStress = EngStress*(1+EngStrain)
    plt.figure(3)
    plt.grid("on")
    plt.title("True Strain VS. True Stress Curve")
    plt.xlabel("True Strain[-]")
    plt.ylabel("True Stress[GPa]")
    plt.plot(TrueStrain,TrueStress)
    plt.legend(["0.01mm/ms","0.1mm/ms","1mm/ms","10mm/ms"])
    return TrueStrain,TrueStress,[PeakStrain,PeakStress]

def MAT24LCSS(x,y,yieldpoint=0.02,Ymould_user=0,Peak=[0,0],strainturn=0,scaleturn=1,curvescale=1,ratio=1,extend=0,scaleextend=0,alignstrain=1.45,pointnum=50,smoothflag=0):
    '''Input:
        TrueStrain,TrueStressm
        Output:Modiifed Eff. Stress Vs. Eff. plastic strain
    '''
    y_strain = []
    x_strain = []
    for index,value in enumerate(x):
        if value >= yieldpoint and (value > x[index-1] or y[index]>y[index-1]):
            x_strain.append(value - y[index]/Ymould_user)
            y_strain.append(y[index])
    plt.figure(4)
    displayratio = int(len(x_strain)*ratio)
    x_strain_mod =[0] 
    y_strain_mod = [y_strain[0]/3.]
    delta = int(len(x_strain)/pointnum)
    x_strain_mod.extend( [value for index,value in enumerate(x_strain) if index%delta==0]);y_strain_mod.extend([value for index,value in enumerate(y_strain) if index%delta==0])
    for index,value in enumerate(x_strain_mod):
        if value > strainturn:
            y_strain_mod[index]= y_strain_mod[index-1]+scaleturn*(x_strain_mod[index]-x_strain_mod[index-1])
        else:
            continue
    if x_strain_mod[-1] < alignstrain:
        x_strain_mod.append(alignstrain)
        y_strain_mod.append(y_strain_mod[-1]+scaleturn*(x_strain_mod[-1]-x_strain_mod[-2]))
    for index,value in enumerate(x_strain_mod):
        y_strain_mod[index] = y_strain_mod[index]*curvescale
    if extend != 0:
        x_strain_mod.append(x_strain_mod[-1]+extend)
        y_strain_mod.append(y_strain_mod[-1]+extend*(scaleextend))
        plt.title("Modified Effective Plastic Strain VS. Stress Curve")
    plt.xlabel("Effective Plastic Strain[-]")
    plt.ylabel("Stress[GPa]")
    plt.grid('on')
    if smoothflag == 1:
        print "The curve has been smoothed!"
        f = interpolate.interp1d(x_strain_mod,y_strain_mod,kind='cubic')
        x_strain_smooth = np.linspace(x_strain_mod[0],x_strain_mod[-1],pointnum)
        y_strain_smooth = f(x_strain_smooth)
        plt.plot(x_strain_smooth[:displayratio],y_strain_smooth[:displayratio])
#        plt.legend(["Strain Rate 0.001/ms","Modified Strain Rate 0.001/ms","Strain Rate 0.01/ms","Modified  Strain Rate 0.01/ms","Strain Rate 0.1/ms","Modified Strain Rate 0.1/ms","Strain Rate 1/ms","Modified Strain Rate 1/ms",])
        return x_strain_smooth,y_strain_smooth
    elif smoothflag == 0:
        plt.plot(x_strain_mod[:displayratio],y_strain_mod[:displayratio])
#        plt.legend(["Strain Rate 0.001/ms","Modified Strain Rate 0.001/ms","Strain Rate 0.01/ms","Modified  Strain Rate 0.01/ms","Strain Rate 0.1/ms","Modified Strain Rate 0.1/ms","Strain Rate 1/ms","Modified Strain Rate 1/ms",])
        return x_strain_mod,y_strain_mod

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
    return x_fit,y_fit,x,y
#    plt.show()

def CowperSymondsCurve(File,A,l0,Ymould_user,pointnum,ratelist,x_p,y_p):
   y_strain = []
   x_strain = []
   for index,value in enumerate(x):
       if value >= yieldpoint and (value > x[index-1] or y[index]>y[index-1]):
           x_strain.append(value - y[index]/Ymould_user)
           y_strain.append(y[index])
   plt.figure(4)
   plt.plot(x_strain,y_strain,'--')
   popt,pcov = optimize.curve_fit(CowperSymondsFunc,x_p,y_p,[0.001,0.001])
   c = popt[0];p = popt[1]
#    print "Cowper-Symonds Constants:\nc:%f\tp:%f\n" %(c,p)
   x_strain_fit = np.linspace(0,x_strain[-1],pointnum)
   for rate in ratelist:
       plt.figure(4)
       y_strain_fit = CowperSymondsFunc(x_strain_fit,c,p,rate=0.001)
       plt.plot(x_strain_fit,y_strain_fit)
   return x_strain_fit,y_strain_fit

def DSGZ_CurveKey(Ymould_user,T,K,c1,c2,c3,c4,a,m,alpha,strainrate_list,curvenum,num):
    res_stress = []
    res_strain = []
    # x,y = CowperSymondsCurve(FDfile[0],A,l0,Ymould_user,80,ratelist,x_p,y_p)
    fout = open("Curve.key",'w')
    fout.write('*KEYWORD\n')
    fout.write( "$%10s%10s%10s%10s%10s%10s%10s%10s\n" %('c1','c2','m','a','K','c3','c4','alpha'))
    fout.write( "$%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n" %(c1,c2,m,a,K,c3,c4,alpha))
    fout.write('$ Created: ' + time.strftime("%d.%m.%Y %H:%M:%S") + '\n')
    fout.write("*DEFINE_TABLE\n%d\n" %(curvenum))
    for i in strainrate_list:
        fout.write("%f\n" %(i))
    for index,strainrate in enumerate(strainrate_list):
        curvenum += 1
        strain = np.linspace(0,2,num)
        fout.write('*DEFINE_CURVE_TITLE\nRate %.5f\n' %(strainrate_list[index]))
        fout.write('$     LCID      SIDR       SFA       SFO      OFFA      OFFO    DATTYP\n')
        fout.write('      %d         0    1.0000    1.0000    0.0000    0.0000\n' %(curvenum))
        print strainrate
        stress = stress_strain_strainrate_T(strain,strainrate,T,K,c1,c2,c3,c4,a,m,alpha)
        res_stress.append(stress)
#        res_strain.append(strain)
        res_strain.append([strain[i] - stress[i]/Ymould_user for i in range(num) ])
        for i in range(num):
            if res_strain[index][i] > 0:
                fout.write("%f,%f\n" %(res_strain[index][i],res_stress[index][i]))
    fout.write("*END\n")
    fout.close()
    return res_strain,res_stress
    
if __name__ == '__main__':
    A = 30; l0 = 10 #样件截面积，标距
    PATH = r'C:\Users\chaos\Desktop\paper\paper\HT\\'
    FDfile = ["test100.txt","test1000.txt","test10000.txt"]#力位移曲线
    FDfile = [PATH + i for i in FDfile]
    strain_rate = [-6.907,-4.605,-2.302,0] #力位移曲线对应的应变率
    log_strain_rate = [np.power(np.e,i) for i in strain_rate]
    curvenum = 2350 #key文件中曲线的编号起始编号
    ratio = 1 #显示比例
#################################User Define###################################
    Ymould_user = 1.5 #定义弹性模量
    yieldpoint= 0.100 #屈服点（避开交叉区域）
    pointnum = 80     #输出点的个数
    strainturnlist=[1.3,1.3,1.35,0.9]#指定每根曲线发生转折的应变
    scaleturnlist=[0.08,0.08,0.06,0.0450];#指定每根曲线发生转折的比例
    curvescalelist = [1,1.05,1.05,1.1]#指定每个区县整体偏移的比例
#    strainturnlist=[1.2,1.2,1.2,0.75];scaleturnlist=[0.06,0.06,0.055,0.050];curvescalelist = [1,1.,1.05,1.05]
    alignstrain = 1.45#应力应变曲线对齐的位置（应变值）
    extend = 0.5     #对齐后曲线延长距离（应变值）
    scaleextend=0.08 #延长段的增量
#################################User Define###################################
    res = CurveKey(A,l0,FDfile,strain_rate,curvenum,ratio,Ymould_user,yieldpoint,pointnum,strainturnlist,scaleturnlist,curvescalelist,alignstrain,extend,scaleextend)

    FILE_list = ['bndout100','bndout1000','bndout10000']
    loadspeed_list = [0.1,1.,10.]
    plt.figure(1)

    for FILE,testfile,loadspeed in zip(FILE_list,["test100.txt","test1000.txt","test10000.txt"],loadspeed_list):
        BndoutFile = alv.bndout(PATH+FILE,0,0,0,0)
        b = BndoutFile.run
        num_point = b['time'].index[-1]
        print '1111',num_point
        try:
            neg_pos = b[b['yforce']<0].index[0] 
            b['yforce'].iloc[neg_pos] = 0
        except:
            neg_pos =  num_point
        x_sim = b['time'][2:neg_pos+1]*loadspeed
        y_sim = b['yforce'][2:neg_pos+1]
        plt.plot(x_sim,y_sim,'--')
        data = pd.read_csv(PATH+testfile)
        x_exp = data.Disp
        y_exp = data.Force
        y_exp =  pd.Series(signal.medfilt(y_exp,201))
        num_exp_point = y_exp.index[-1]
        print '2222',num_exp_point
        delta_exp = int(float(num_exp_point)/num_point)
        y_exp_cor = [i for i in y_exp if i%delta_exp==0]
        print len(list(y_exp_cor))
        print np.trapz(b['yforce'][2:neg_pos+1],b['time'][2:neg_pos+1]*loadspeed)
        print np.corrcoef(y_exp_cor,y_sim)
    # plt.plot(bb,b['yforce'])
    plt.legend(["0.1mm/ms exp.","1mm/ms exp.","10mm/ms exp.","0.1mm/ms prediction","1mm/ms prediction","10mm/ms prediction"],loc='lower right')
    plt.grid('on')
    plt.show()
    
    

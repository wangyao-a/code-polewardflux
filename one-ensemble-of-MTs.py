#!/usr/bin/env python
# coding: utf-8



import random
import math
import numpy as np
import pandas as pd

def stepnum(ypsl0,ypsl1,ypsl2):#Determine if each pair of MT-bound heads of kinesin-5 motor step forward, step backward
    kpos = 12.8
    kneg = kpos/15
    ED=2.5
    kNL=200
    kD=50
    kr=4.6
    p=1
    P0T=(1/kD+1/kr)/(1/kD+1/kr+1/kpos)
    P0L=(1/kD+1/kr)/(1/kD+1/kr+1/kneg)
    PE1=math.exp(ED+((ypsl0-ypsl1)*p)/4.11)/(math.exp(ED+((ypsl0-ypsl1)*p)/4.11)+1)
    PE11=1-PE1
    PE2=math.exp(ED+((ypsl2-ypsl0)*p)/4.11)/(math.exp(ED+((ypsl2-ypsl0)*p)/4.11)+1)
    PE21=1-PE2
    kT=PE1/(1/kpos+P0T/kr)+PE11/(1/kD+1/kpos+1/kr)
    kL=PE2/(1/kD+1/kr+1/kneg)+PE21/(1/kneg+P0L/kr)
    ran3 = random.uniform(0, 1)
    ran4 = random.uniform(0, 1)
    na1 = 0
    if ran3 < kT*PE1 * h:
        na1+= 1
        #print(kT*PE1,kL*PE21,ypsl0,ypsl1,ypsl2)
    if ran4 < kL*PE21 * h:
        na1-= 1
        #print(kT*PE1,kL*PE21,ypsl0,ypsl1,ypsl2)
    return na1


# In[2]:


def rebind():#Determine if another pair of detached heads of the MT-bound kinesin-5 bind to MT.
    u=0.2
    ran4 = random.uniform(0, 1)
    Kin=0
    if ran4 < u * h:
        Kin = 1
    return Kin        


# In[3]:


def dissociation(F,ypsl0,ypsl1,ypsl2):#Determine if each pair of MT-bound heads of kinesin-5 motor detach. 
    F=-F
    kpos = 12.8
    kneg = kpos/15
    ED=2.5
    kNL=200
    kD=50
    dpos=8.2
    delta=1
    kr=4.6
    kw0=5
    p=1
    es0=0.1
    deltas=2.3
    ran3 = random.uniform(0, 1)
    if F <= 0 :
        P1=kneg/(kNL+kneg)
    if F >= 4:
        P1=kpos/(kNL+kpos)
    if F > 0 and F < 4:
        k=kneg+(kpos-kneg)*F/4
        P1=k/(kNL+k)            
    Pd1=1
    PE1=math.exp(ED+((ypsl0-ypsl1)*p)/4.11)/(math.exp(ED+((ypsl0-ypsl1)*p)/4.11)+1)
    PE11=1-PE1
    PE2=math.exp(ED+((ypsl2-ypsl0)*p)/4.11)/(math.exp(ED+((ypsl2-ypsl0)*p)/4.11)+1)
    PE21=1-PE2
    P2T=PE1*kpos/(kpos+kD)+PE11*kneg/(kneg+kD)
    P2L=PE2*kpos/(kpos+kD)+PE21*kneg/(kneg+kD)
    kd2=kw0*math.exp(abs(F)*delta/4.130238407)
    Pd2=kd2/(kd2+kD)
    P0T=(1/kD+1/kr)/(1/kD+1/kr+1/kpos)
    P0L=(1/kD+1/kr)/(1/kD+1/kr+1/kneg)
    kT=PE1/(1/kpos+P0T/kr)+(1-PE1)/(1/kD+1/kpos+1/kr)
    kL=PE2/(1/kD+1/kr+1/kneg)+(1-PE2)/(1/kneg+P0L/kr)
    episilon=kT*P1*Pd1+(kT*P2T+kL*P2L)*Pd2+es0*math.exp(abs(F)*deltas/4.130238407)
    Kin=1
    #print(episilon)
    if ran3 < episilon * h:
        Kin= 0
    return Kin


# In[4]:


def jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx):#calculate the displacements of MTs, spindle poles, and kinetochores to equilibrate the spindle apparatus again
    mtui=[0 for i in range(8)]
    num1=0
    num2=0 
    nump1=0
    nump2=0
    numpx1=0
    numpx2=0
    F0=Kp1*(l1l-Xpole)+Kp2*(l2l-Xpole)
    F7=Kp1*(l4r-Xpo)+Kp2*(l3r-Xpo)
    F5=Kp3*(l1r-Xkin)+Kp4*(Xkin1-Xkin-kdis)
    F6=Kp4*(Xkin+kdis-Xkin1)+Kp3*(l4l-Xkin1)
    if abs(F0)<0.0001 and abs(F1)<0.00001 and abs(F2)<0.00001 and abs(F3)<0.00001 and abs(F4)<0.00001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
        return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]
    
    nu=len(kina)
    nu1=len(kinap)
    nu1x=len(kinapx)
    for ton1 in range(nu):
        if kina[ton1]==1 and kinb[ton1]==1:
            if dxc[ton1]>xdmotor:
                num1+=1
            if dxc[ton1]<-xdmotor:
                num2+=1
    for ton1 in range(nu1):
        if kinap[ton1]==1 and kinbp[ton1]==1:
            if dxcp[ton1]>xdmotor:
                nump1+=1
            if dxcp[ton1]<-xdmotor:
                nump2+=1
    for ton1 in range(nu1x):
        if kinapx[ton1]==1 and kinbpx[ton1]==1:
            if dxcpx[ton1]>xdmotor:
                numpx1+=1
            if dxcpx[ton1]<-xdmotor:
                numpx2+=1
    jisuan=sorted(dxc)
    for ton1 in range(nu-monum):
        jisuan.remove(0)
    jisuanp=sorted(dxcp)
    for ton1 in range(nu1-monump):
        jisuanp.remove(0)
    jisuanpx=sorted(dxcpx)
    for ton1 in range(nu1x-monumpx):
        jisuanpx.remove(0)
    cnn0=0
    zq=0
    bd=monum-num1-num2
    bdp=monump-nump1-nump2
    bdpx=monumpx-numpx1-numpx2

    m1=[monum-num1]
    m0=[num2]
    m2=[0]
    if monum>0:
        ttt=jisuan[0]
        for tt1 in range(monum):
            if abs(jisuan[tt1]-ttt)>0.0001:
                if jisuan[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuan[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuan[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuan[tt1]
        m2.append(num2)
        m0.append(monum-num1)
        m1.append(monum)
    fen2=list(set(m2))
    fen1=list(set(m1))
    fen0=list(set(m0))
    fen2.sort(key=m2.index)
    fen1.sort(key=m1.index)
    fen0.sort(key=m0.index)

    m1=[monump-nump1]
    m0=[nump2]
    m2=[0]
    if monump>0:
        ttt=jisuanp[0]
        for tt1 in range(monump):
            if abs(jisuanp[tt1]-ttt)>0.0001:
                if jisuanp[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuanp[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuanp[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuanp[tt1]
        m2.append(nump2)
        m0.append(monump-nump1)
        m1.append(monump)
    fenp2=list(set(m2))
    fenp1=list(set(m1))
    fenp0=list(set(m0))
    fenp2.sort(key=m2.index)
    fenp1.sort(key=m1.index)
    fenp0.sort(key=m0.index)

    m1=[monumpx-numpx1]
    m0=[numpx2]
    m2=[0]
    if monumpx>0:
        ttt=jisuanpx[0]
        for tt1 in range(monumpx):
            if abs(jisuanpx[tt1]-ttt)>0.0001:
                if jisuanpx[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuanpx[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuanpx[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuanpx[tt1]
        m2.append(numpx2)
        m0.append(monumpx-numpx1)
        m1.append(monumpx)
    fenpx2=list(set(m2))
    fenpx1=list(set(m1))
    fenpx0=list(set(m0))
    fenpx2.sort(key=m2.index)
    fenpx1.sort(key=m1.index)
    fenpx0.sort(key=m0.index)
    ########################################
    dzk=[[0 for xxn2 in range(len(fen0))] for xxn1 in range(2)]
    dz0=[[0 for xxn2 in range(len(fen2))],[0 for xxn1 in range(len(fen1))]]
    dzkp=[[0 for xxn2 in range(len(fenp0))] for xxn1 in range(2)]
    dz0p=[[0 for xxn2 in range(len(fenp1))],[0 for xxn1 in range(len(fenp2))]]
    dzkpx=[[0 for xxn2 in range(len(fenpx0))] for xxn1 in range(2)]
    dz0px=[[0 for xxn2 in range(len(fenpx1))],[0 for xxn1 in range(len(fenpx2))]]
    ta1=-1
    for xck in fen0:
        ta1+=1
        cnnk=monum-num1-fen0[-ta1-1]
        chek=monum-num1-1
        for xxxn in range(cnnk):
            dzk[0][ta1]+=jisuan[chek-xxxn]-xdmotor
    ta1=-1
    for xck in fen0:
        ta1+=1
        cnnk=xck-num2
        chek=num2
        for xxxn in range(cnnk):
            dzk[1][ta1]+=jisuan[chek+xxxn]+xdmotor 
    ta2=-1
    for xckp in fenp0:
        ta2+=1
        cnnkp=xckp-nump2
        chekp=nump2
        for xxxn in range(cnnkp):
            dzkp[0][ta2]+=jisuanp[chekp+xxxn]+xdmotor
    ta2=-1
    for xckp in fenp0:
        ta2+=1
        cnnkp=monump-nump1-fenp0[-ta2-1]
        chekp=monump-nump1-1 
        for xxxn in range(cnnkp):
            dzkp[1][ta2]+=jisuanp[chekp-xxxn]-xdmotor 
    ta2=-1
    for xckpx in fenpx0:
        ta2+=1
        cnnkpx=xckpx-numpx2
        chekpx=numpx2
        for xxxn in range(cnnkpx):
            dzkpx[0][ta2]+=jisuanpx[chekpx+xxxn]+xdmotor
    ta2=-1
    for xckpx in fenpx0:
        ta2+=1
        cnnkpx=monumpx-numpx1-fenpx0[-ta2-1]
        chekpx=monumpx-numpx1-1 
        for xxxn in range(cnnkpx):
            dzkpx[1][ta2]+=jisuanpx[chekpx-xxxn]-xdmotor 
    ta3=-1
    for xc0 in fen2:
        ta3+=1
        cnn0=num2-fen2[-ta3-1]
        che0=num2-1
        for xxxn in range(cnn0):
            dz0[0][ta3]+=jisuan[che0-xxxn]+xdmotor
    ta3=-1
    for xc0 in fen1:    
        ta3+=1
        cnn0=xc0-monum+num1
        che0=monum-num1
        for xxxn in range(cnn0):
            dz0[1][ta3]+=jisuan[che0+xxxn]-xdmotor
    ta4=-1
    for xc0p in fenp1:
        ta4+=1
        cnn0p=xc0p-monump+nump1
        che0p=monump-nump1
        for xxxn in range(cnn0p):
            dz0p[0][ta4]+=jisuanp[che0p+xxxn]-xdmotor
    ta4=-1
    for xc0p in fenp2:
        ta4+=1
        cnn0p=nump2-fenp2[-ta4-1]
        che0p=nump2-1
        for xxxn in range(cnn0p):
            dz0p[1][ta4]+=jisuanp[che0p-xxxn]+xdmotor
    ta4=-1
    for xc0px in fenpx1:
        ta4+=1
        cnn0px=xc0px-monumpx+numpx1
        che0px=monumpx-numpx1
        for xxxn in range(cnn0px):
            dz0px[0][ta4]+=jisuanpx[che0px+xxxn]-xdmotor
    ta4=-1
    for xc0px in fenpx2:
        ta4+=1
        cnn0px=numpx2-fenpx2[-ta4-1]
        che0px=numpx2-1
        for xxxn in range(cnn0px):
            dz0px[1][ta4]+=jisuanpx[che0px-xxxn]+xdmotor
    a00=-Kp1-Kp2;a10=Kp2;a30=Kp1
    a21=Kp2;a41=Kp1;a71=-Kp1-Kp2
    a32=Kp3;a52=-Kp3-Kp4;a62=Kp4
    a43=Kp3;a53=Kp4;a63=-Kp3-Kp4 
   
    for xx1 in range(2):
        ta1=-1
        for xck in fen0:
            ta1+=1
            if xx1==0: 
                if ta1<len(dzk[0]):
                    cnnk=monum-num1-fen0[-ta1-1]
                if ta1>=len(dzk[0]):
                    break
            if xx1==1:
                if ta1<len(dzk[1]):
                    cnnk=xck-num2
                if ta1>=len(dzk[1]):
                    break
            for xx2 in range(2):
                ta2=-1
                for xckp in fenp0:
                    ta2+=1
                    if xx2==0:
                        if ta2<len(dzkp[0]):
                            cnnkp=xckp-nump2
                        if ta2>=len(dzkp[0]):
                            break
                    if xx2==1:
                        if ta2<len(dzkp[1]):
                            cnnkp=monump-nump1-fenp0[-ta2-1]
                        if ta2>=len(dzkp[1]):
                            break
                    for xx3 in range(2):
                        ta5=-1
                        for xckpx in fenpx0:    
                            ta5+=1
                            if xx3==0:
                                if ta5<len(dzkpx[0]):
                                    cnnkpx=xckpx-numpx2
                                if ta5>=len(dzkpx[0]):
                                    break
                            if xx3==1:
                                if ta5<len(dzkpx[1]):
                                    cnnkpx=monumpx-numpx1-fenpx0[-ta5-1]
                                if ta5>=len(dzkpx[1]):
                                    break
                            if xx1==0:
                                xa3=fen2
                            if xx1==1:
                                xa3=fen1
                            ta3=-1
                            for xc0 in xa3:
                                #print(xqq,xc0,len(dz0[xqq][1]),len(dz0[xqq][0]),'3')
                                ta3+=1  
                                if xx1==0:
                                    if ta3<len(dz0[0]):
                                        cnn0=num2-fen2[-ta3-1]
                                    if ta3>=len(dz0[0]):
                                        break
                                if xx1==1:
                                    if ta3<len(dz0[1]):
                                        cnn0=xc0-monum+num1
                                    if ta3>=len(dz0[1]):
                                        break
                                if xx2==0:
                                    xa4=fenp1
                                if xx2==1:
                                    xa4=fenp2
                                ta4=-1
                                for xc0p in xa4:
                                    ta4+=1
                                    if xx2==0:
                                        if ta4<len(dz0p[0]):
                                            cnn0p=xc0p-monump+nump1
                                        if ta4>=len(dz0p[0]):
                                            break
                                    if xx2==1:
                                        if ta4<len(dz0p[1]):
                                            cnn0p=nump2-fenp2[-ta4-1]    
                                    if xx3==0:
                                        xa5=fenpx1
                                    if xx3==1:
                                        xa5=fenpx2
                                    ta6=-1
                                    for xc0px in xa5:
                                        ta6+=1
                                        if xx3==0:
                                            if ta6<len(dz0px[0]):
                                                cnn0px=xc0px-monumpx+numpx1
                                            if ta6>=len(dz0px[0]):
                                                break
                                        if xx3==1:
                                            if ta6<len(dz0px[1]):
                                                cnn0px=numpx2-fenpx2[-ta6-1] 
                                            if ta6>len(dz0px[1]):
                                                break
                                        a04=-Kp1;a24=round(-(nump1+nump2+cnnkp-cnn0p)*Km,3);a34=round((nump1+nump2+cnnkp-cnn0p)*Km+Kp1+Kp3,3);a54=-Kp3
                                        a05=-Kp2;a15=round(((numpx1+numpx2+cnnkpx-cnn0px)+(num1+num2+cnnk-cnn0))*Km+Kp2,3);
                                        a25=round(-Km*(num1+num2+cnnk-cnn0),3);a45=round(-(numpx1+numpx2+cnnkpx-cnn0px)*Km,3)
                                        a16=round((num1+num2+cnnk-cnn0)*Km,3);a26=round(-Km*(nump1+nump2+cnnkp-cnn0p+num1+num2+cnnk-cnn0)-Kp2,3)
                                        a36=round(Km*(nump1+nump2+cnnkp-cnn0p),3);a76=Kp2
                                        a17=round(Km*(numpx1+numpx2+cnnkpx-cnn0px),3);a47=round(-Km*(numpx1+numpx2+cnnkpx-cnn0px)-Kp1-Kp3,3);a67=Kp3;a77=Kp1
                                        b4=-Km*(dzkp[xx2][ta2]-dz0p[xx2][ta4])-F1
                                        b5=-Km*(-dzkpx[xx3][ta5]+dz0px[xx3][ta6]+dzk[xx1][ta1]-dz0[xx1][ta3])-F2
                                        b6=-Km*(dzkp[xx2][ta2]-dz0p[xx2][ta4]+dzk[xx1][ta1]-dz0[xx1][ta3])-F3
                                        b7=Km*(dzkpx[xx3][ta5]-dz0px[xx3][ta6])-F4
                                        b=np.array([round(-F0,14),round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14)])
                                            
                                        try:
                                            A=np.array([[a00,a10,0,a30,0,0,0,0],[0,0,a21,0,a41,0,0,a71],[0,0,0,a32,0,a52,a62,0],
                                                        [0,0,0,0,a43,a53,a63,0],[a04,0,a24,a34,0,a54,0,0],[a05,a15,a25,0,a45,0,0,0],
                                                        [0,a16,a26,a36,0,0,0,a76],[0,a17,0,0,a47,0,a67,a77]])
                                            #b=np.array([-F0,-F7,-F5,-F6,b4,b5,b6,b7])

                                            mtui=np.linalg.solve(A,b)
                                        except:
                                            A=np.array([[0,a21,0,a41,0,0,a71],[0,0,a32,0,a52,a62,0],
                                                    [0,0,0,a43,a53,a63,0],[0,a24,a34,0,a54,0,0],[a15,a25,0,a45,0,0,0],
                                                    [a16,a26,a36,0,0,0,a76],[a17,0,0,a47,0,a67,a77]])
                                            b=np.array([round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14)])
                                          
                                            xtui=np.linalg.solve(A,b)
                                            for xn in range(7):
                                                mtui[xn+1]=xtui[xn]
                                            mtui[0]=0
                                        tuzx=(mtui[0]+mtui[7])/2
                                        for xn in range(8):
                                            mtui[xn]-=tuzx

                                        check0=-1#################
                                        checkk=-1################
                                        checkp0=-1#################
                                        checkpk=-1#################
                                        checkpx0=-1#################
                                        checkpxk=-1#################
                                        if xx1==0 and mtui[1]-mtui[2]>=0:#
                                            check0=0#################
                                            checkk=0#################
                                            for ton1 in range(monum):
                                                if jisuan[ton1]<-xdmotor and abs(jisuan[ton1]+mtui[1]-mtui[2])<=xdmotor:
                                                    check0+=1
                                                if abs(jisuan[ton1])<=xdmotor and jisuan[ton1]+mtui[1]-mtui[2]>xdmotor:
                                                    checkk+=1
                                        if xx1==1 and mtui[1]-mtui[2]<=0:
                                            check0=0#################
                                            checkk=0#################
                                            for ton1 in range(monum):
                                                if abs(jisuan[ton1])<=xdmotor and jisuan[ton1]+mtui[1]-mtui[2]<-xdmotor:
                                                    checkk+=1
                                                if jisuan[ton1]>xdmotor and abs(jisuan[ton1]+mtui[1]-mtui[2])<=xdmotor:
                                                    check0+=1
                                        if check0!=cnn0 or checkk!=cnnk:
                                            continue
                                        if xx2==1 and mtui[3]-mtui[2]>=0:
                                            checkp0=0#################
                                            checkpk=0#################
                                            for ton1 in range(monump):
                                                if jisuanp[ton1]<-xdmotor and abs(jisuanp[ton1]+mtui[3]-mtui[2])<=xdmotor:
                                                    checkp0+=1
                                                if abs(jisuanp[ton1])<=xdmotor and jisuanp[ton1]+mtui[3]-mtui[2]>xdmotor:
                                                    checkpk+=1
                                        if xx2==0 and mtui[3]-mtui[2]<=0:#
                                            checkp0=0#################
                                            checkpk=0#################
                                            for ton1 in range(monump):
                                                if abs(jisuanp[ton1])<=xdmotor and jisuanp[ton1]+mtui[3]-mtui[2]<-xdmotor:
                                                    checkpk+=1
                                                if jisuanp[ton1]>xdmotor and abs(jisuanp[ton1]+mtui[3]-mtui[2])<=xdmotor:
                                                    checkp0+=1
                                        if checkp0!=cnn0p or checkpk!=cnnkp:
                                            continue      
                                        if xx3==1 and mtui[4]-mtui[1]>=0:
                                            checkpx0=0#################
                                            checkpxk=0#################
                                            for ton1 in range(monumpx):
                                                if jisuanpx[ton1]<-xdmotor and abs(jisuanpx[ton1]+mtui[4]-mtui[1])<=xdmotor:
                                                    checkpx0+=1
                                                if abs(jisuanpx[ton1])<=xdmotor and jisuanpx[ton1]+mtui[4]-mtui[1]>xdmotor:
                                                    checkpxk+=1
                                        if xx3==0 and mtui[4]-mtui[1]<=0:#
                                            checkpx0=0#################
                                            checkpxk=0#################
                                            for ton1 in range(monumpx):
                                                if abs(jisuanpx[ton1])<=xdmotor and jisuanpx[ton1]+mtui[4]-mtui[1]<-xdmotor:
                                                    checkpxk+=1
                                                if jisuanpx[ton1]>xdmotor and abs(jisuanpx[ton1]+mtui[4]-mtui[1])<=xdmotor:
                                                    checkpx0+=1
                                        if checkpx0!=cnn0px or checkpxk!=cnnkpx:
                                            continue    
                                        #print(mtui,'qqq')
                                        return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]
    tq=1
    if abs(F2)>0.0001 and abs(F3)<0.0001:
        tq=0
    lastui1=-1000
    if tq==0:
        for xn in range(num1+numpx2+num2+numpx1+1):
            if xn < num1:
                tui1=-(jisuan[-num1+xn]-xdmotor+0.1)
            if xn >= num1 and xn <num1+numpx2:
                tui1=(jisuanpx[numpx2-1-xn+num1]+xdmotor-0.1)
            if xn >= num1+numpx2 and xn <num1+numpx2+num2:
                tui1=-(jisuan[num2-1-xn+num1+numpx2]+xdmotor-0.1)
            if xn >= num1+numpx2+num2 and xn <num1+numpx2+num2+numpx1:
                tui1=(jisuanpx[monumpx-numpx1+xn-num1-numpx2-num2]-xdmotor+0.1)
            if xn==num1+numpx2+num2+numpx1:
                tui1=0
            if abs(tui1-lastui1)<0.0001:
                continue
            lastui1=tui1
            lastui2=-1000
            for xn1 in range(nump2+nump1+1):#MT4
                if xn1 <nump1:
                    tui2=(-jisuanp[-nump1+xn1]+xdmotor-0.1)
                if xn1 >= nump1 and xn1 <nump1+nump2:
                    tui2=-(jisuanp[nump2-1-xn1+nump1]+xdmotor-0.1)
                if xn1==nump1+nump2:
                    tui2=0
                if abs(tui2-lastui2)<0.0001 or abs(tui1)+abs(tui2)==0:
                    continue
                lastui2=tui2
                #print(tui1,tui2)
                F0=Kp1*(l1l+tui2-Xpole)+Kp2*(l2l+tui1-Xpole)
                F7=Kp1*(l4r-Xpo)+Kp2*(l3r-Xpo)
                F5=Kp3*(l1r+tui2-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                F6=Kp4*(Xkin+kdis-Xkin1)+Kp3*(l4l-Xkin1)
                F3=Kp2*(Xpo-l3r)
                F1=Kp1*(l1l+tui2-Xpole)+Kp3*(l1r+tui2-Xkin)
                F2=Kp2*(l2l+tui1-Xpole)
                F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                for ton1 in range(monum):
                    if jisuan[ton1]+tui1>xdmotor:
                        F3+=(jisuan[ton1]+tui1-xdmotor)*Km
                        F2+=(jisuan[ton1]+tui1-xdmotor)*Km
                    if jisuan[ton1]+tui1<-xdmotor:
                        F3+=(jisuan[ton1]+tui1+xdmotor)*Km
                        F2+=(jisuan[ton1]+tui1+xdmotor)*Km
                for ton1 in range(monump):
                    if jisuanp[ton1]+tui2>xdmotor:
                        F3+=(jisuanp[ton1]+tui2-xdmotor)*Km
                        F1+=(jisuanp[ton1]+tui2-xdmotor)*Km
                    if jisuanp[ton1]+tui2<-xdmotor:
                        F3+=(jisuanp[ton1]+tui2+xdmotor)*Km
                        F1+=(jisuanp[ton1]+tui2+xdmotor)*Km

                for ton1 in range(monumpx):
                    if jisuanpx[ton1]-tui1>xdmotor:
                        F4-=(jisuanpx[ton1]-tui1-xdmotor)*Km
                        F2-=(jisuanpx[ton1]-tui1-xdmotor)*Km
                    if jisuanpx[ton1]-tui1<-xdmotor:
                        F4-=(jisuanpx[ton1]-tui1+xdmotor)*Km
                        F2-=(jisuanpx[ton1]-tui1+xdmotor)*Km
                if abs(F0)<0.0001 and abs(F1)<0.00001 and abs(F2)<0.00001 and abs(F3)<0.00001 and abs(F4)<0.00001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
                    mtui[1]+=tui1
                    mtui[3]+=tui2
                    return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]

                jisuan1=[]
                jisuanp1=[]
                jisuanpx1=[]
                for ton1 in range(monum):
                    jisuan1.append(jisuan[ton1]+tui1)
                for ton1 in range(monump):
                    jisuanp1.append(jisuanp[ton1]+tui2)
                for ton1 in range(monumpx):
                    jisuanpx1.append(jisuanpx[ton1]-tui1)
                tuf,mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]=jiaocheck(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7)
                if tuf==1:
                    mtui[1]+=tui1
                    mtui[3]+=tui2
                    return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]
                    break
    if tq==1:
        for xn in range(num1+nump2+num2+nump1+1):
            if xn < num1:
                tui1=(jisuan[-num1+xn]-xdmotor+0.1)
            if xn >= num1 and xn <num1+nump2:
                tui1=(jisuanp[nump2-1-xn+num1]+xdmotor-0.1)
            if xn >= num1+nump2 and xn <num1+nump2+num2:
                tui1=(jisuan[num2-1-xn+num1+nump2]+xdmotor-0.1)
            if xn >= num1+nump2+num2 and xn <num1+nump2+num2+nump1:
                tui1=(jisuanp[monump-nump1+xn-num1-nump2-num2]-xdmotor+0.1)
            if xn==num1+nump2+num2+nump1:
                tui1=0
            if abs(tui1-lastui1)<0.0001:
                continue
            lastui1=tui1
            lastui2=-1000
            for xn1 in range(numpx2+numpx1+1):#MT4
                if xn1 <numpx1:
                    tui2=(-jisuanpx[-numpx1+xn1]+xdmotor-0.1)
                if xn1 >= numpx1 and xn1 <numpx1+numpx2:
                    tui2=-(jisuanpx[numpx2-1-xn1+numpx1]+xdmotor-0.1)
                if xn1==numpx1+numpx2:
                    tui2=0
                if abs(tui2-lastui2)<0.0001 or abs(tui1)+abs(tui2)==0:
                    continue
                lastui2=tui2
                #print(tui1,tui2)
                F0=Kp1*(l1l-Xpole)+Kp2*(l2l-Xpole)
                F7=Kp1*(l4r+tui2-Xpo)+Kp2*(l3r+tui1-Xpo)
                F5=Kp3*(l1r-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                F6=Kp4*(Xkin+kdis-Xkin1)+Kp3*(l4l+tui2-Xkin1)
                F3=Kp2*(Xpo-l3r-tui1)
                F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                F2=Kp2*(l2l-Xpole)
                F4=Kp1*(Xpo-l4r-tui2)+Kp3*(Xkin1-l4l-tui2)
                for ton1 in range(monum):
                    if jisuan[ton1]-tui1>xdmotor:
                        F3+=(jisuan[ton1]-tui1-xdmotor)*Km
                        F2+=(jisuan[ton1]-tui1-xdmotor)*Km
                    if jisuan[ton1]-tui1<-xdmotor:
                        F3+=(jisuan[ton1]-tui1+xdmotor)*Km
                        F2+=(jisuan[ton1]-tui1+xdmotor)*Km
                for ton1 in range(monump):
                    if jisuanp[ton1]-tui1>xdmotor:
                        F3+=(jisuanp[ton1]-tui1-xdmotor)*Km
                        F1+=(jisuanp[ton1]-tui1-xdmotor)*Km
                    if jisuanp[ton1]-tui1<-xdmotor:
                        F3+=(jisuanp[ton1]-tui1+xdmotor)*Km
                        F1+=(jisuanp[ton1]-tui1+xdmotor)*Km

                for ton1 in range(monumpx):
                    if jisuanpx[ton1]+tui2>xdmotor:
                        F4-=(jisuanpx[ton1]+tui2-xdmotor)*Km
                        F2-=(jisuanpx[ton1]+tui2-xdmotor)*Km
                    if jisuanpx[ton1]+tui2<-xdmotor:
                        F4-=(jisuanpx[ton1]+tui2+xdmotor)*Km
                        F2-=(jisuanpx[ton1]+tui2+xdmotor)*Km
                if abs(F0)<0.0001 and abs(F1)<0.00001 and abs(F2)<0.00001 and abs(F3)<0.00001 and abs(F4)<0.00001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
                    mtui[2]+=tui1
                    mtui[4]+=tui2
                    return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]

                jisuan1=[]
                jisuanp1=[]
                jisuanpx1=[]
                for ton1 in range(monum):
                    jisuan1.append(jisuan[ton1]-tui1)
                for ton1 in range(monump):
                    jisuanp1.append(jisuanp[ton1]-tui1)
                for ton1 in range(monumpx):
                    jisuanpx1.append(jisuanpx[ton1]+tui2)
                tuf,mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]=jiaocheck(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7)
                if tuf==1:
                    mtui[2]+=tui1
                    mtui[4]+=tui2
                    return mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]
                    break
 #######################################################################tui           


# In[5]:


def jiaocheck(jisuan,jisuanp,jisuanpx,F0,F1,F2,F3,F4,F5,F6,F7):#Auxiliary Calculation of Displacement
    num3=0
    num4=0 
    nump3=0
    nump4=0
    numpx3=0
    numpx4=0
    for ton1 in range(monum):
        if jisuan[ton1]>xdmotor:
            num3+=1
        if jisuan[ton1]<-xdmotor:
            num4+=1
    for ton1 in range(monump):
        if jisuanp[ton1]>xdmotor:
            nump3+=1
        if jisuanp[ton1]<-xdmotor:
            nump4+=1
    for ton1 in range(monumpx):
        if jisuanpx[ton1]>xdmotor:
            numpx3+=1
        if jisuanpx[ton1]<-xdmotor:
            numpx4+=1
    cnn0=0
    zq=0
    bd=monum-num3-num4
    bdp=monump-nump3-nump4
    bdpx=monumpx-numpx3-numpx4

    m1=[monum-num3]
    m0=[num4]
    m2=[0]
    if monum>0:
        ttt=jisuan[0]
        for tt1 in range(monum):
            if abs(jisuan[tt1]-ttt)>0.0001:
                if jisuan[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuan[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuan[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuan[tt1]
        m2.append(num4)
        m0.append(monum-num3)
        m1.append(monum)
    fen2=list(set(m2))
    fen1=list(set(m1))
    fen0=list(set(m0))
    fen2.sort(key=m2.index)
    fen1.sort(key=m1.index)
    fen0.sort(key=m0.index)

    m1=[monump-nump3]
    m0=[nump4]
    m2=[0]
    if monump>0:
        ttt=jisuanp[0]
        for tt1 in range(monump):
            if abs(jisuanp[tt1]-ttt)>0.0001:
                if jisuanp[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuanp[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuanp[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuanp[tt1]
        m2.append(nump4)
        m0.append(monump-nump3)
        m1.append(monump)
    fenp2=list(set(m2))
    fenp1=list(set(m1))
    fenp0=list(set(m0))
    fenp2.sort(key=m2.index)
    fenp1.sort(key=m1.index)
    fenp0.sort(key=m0.index)

    m1=[monumpx-numpx3]
    m0=[numpx4]
    m2=[0]
    if monumpx>0:
        ttt=jisuanpx[0]
        for tt1 in range(monumpx):
            if abs(jisuanpx[tt1]-ttt)>0.0001:
                if jisuanpx[tt1]<-xdmotor:
                    m2.append(tt1)
                if jisuanpx[tt1]>=-xdmotor and ttt<xdmotor:
                    m0.append(tt1)
                if jisuanpx[tt1]>xdmotor:
                    m1.append(tt1)
                ttt=jisuanpx[tt1]
        m2.append(numpx4)
        m0.append(monumpx-numpx3)
        m1.append(monumpx)
    fenpx2=list(set(m2))
    fenpx1=list(set(m1))
    fenpx0=list(set(m0))
    fenpx2.sort(key=m2.index)
    fenpx1.sort(key=m1.index)
    fenpx0.sort(key=m0.index)
    #########################################
    dzk=[[0 for xxn2 in range(len(fen0))] for xxn1 in range(2)]
    dz0=[[0 for xxn2 in range(len(fen2))],[0 for xxn1 in range(len(fen1))]]
    dzkp=[[0 for xxn2 in range(len(fenp0))] for xxn1 in range(2)]
    dz0p=[[0 for xxn2 in range(len(fenp1))],[0 for xxn1 in range(len(fenp2))]]
    dzkpx=[[0 for xxn2 in range(len(fenpx0))] for xxn1 in range(2)]
    dz0px=[[0 for xxn2 in range(len(fenpx1))],[0 for xxn1 in range(len(fenpx2))]]
    ta1=-1
    for xck in fen0:
        ta1+=1
        cnnk=monum-num3-fen0[-ta1-1]
        chek=monum-num3-1
        for xxxn in range(cnnk):
            dzk[0][ta1]+=jisuan[chek-xxxn]-xdmotor
    ta1=-1
    for xck in fen0:
        ta1+=1
        cnnk=xck-num4
        chek=num4
        for xxxn in range(cnnk):
            dzk[1][ta1]+=jisuan[chek+xxxn]+xdmotor 
    ta2=-1
    for xckp in fenp0:
        ta2+=1
        cnnkp=xckp-nump4
        chekp=nump4
        for xxxn in range(cnnkp):
            dzkp[0][ta2]+=jisuanp[chekp+xxxn]+xdmotor
    ta2=-1
    for xckp in fenp0:
        ta2+=1
        cnnkp=monump-nump3-fenp0[-ta2-1]
        chekp=monump-nump3-1 
        for xxxn in range(cnnkp):
            dzkp[1][ta2]+=jisuanp[chekp-xxxn]-xdmotor 
    ta2=-1
    for xckpx in fenpx0:
        ta2+=1
        cnnkpx=xckpx-numpx4
        chekpx=numpx4
        for xxxn in range(cnnkpx):
            dzkpx[0][ta2]+=jisuanpx[chekpx+xxxn]+xdmotor
    ta2=-1
    for xckpx in fenpx0:
        ta2+=1
        cnnkpx=monumpx-numpx3-fenpx0[-ta2-1]
        chekpx=monumpx-numpx3-1 
        for xxxn in range(cnnkpx):
            dzkpx[1][ta2]+=jisuanpx[chekpx-xxxn]-xdmotor 
    ta3=-1
    for xc0 in fen2:
        ta3+=1
        cnn0=num4-fen2[-ta3-1]
        che0=num4-1
        for xxxn in range(cnn0):
            dz0[0][ta3]+=jisuan[che0-xxxn]+xdmotor
    ta3=-1
    for xc0 in fen1:    
        ta3+=1
        cnn0=xc0-monum+num3
        che0=monum-num3
        for xxxn in range(cnn0):
            dz0[1][ta3]+=jisuan[che0+xxxn]-xdmotor
    ta4=-1
    for xc0p in fenp1:
        ta4+=1
        cnn0p=xc0p-monump+nump3
        che0p=monump-nump3
        for xxxn in range(cnn0p):
            dz0p[0][ta4]+=jisuanp[che0p+xxxn]-xdmotor
    ta4=-1
    for xc0p in fenp2:
        ta4+=1
        cnn0p=nump4-fenp2[-ta4-1]
        che0p=nump4-1
        for xxxn in range(cnn0p):
            dz0p[1][ta4]+=jisuanp[che0p-xxxn]+xdmotor
    ta4=-1
    for xc0px in fenpx1:
        ta4+=1
        cnn0px=xc0px-monumpx+numpx3
        che0px=monumpx-numpx3
        for xxxn in range(cnn0px):
            dz0px[0][ta4]+=jisuanpx[che0px+xxxn]-xdmotor
    ta4=-1
    for xc0px in fenpx2:
        ta4+=1
        cnn0px=numpx4-fenpx2[-ta4-1]
        che0px=numpx4-1
        for xxxn in range(cnn0px):
            dz0px[1][ta4]+=jisuanpx[che0px-xxxn]+xdmotor
    a00=-Kp1-Kp2;a10=Kp2;a30=Kp1
    a21=Kp2;a41=Kp1;a71=-Kp1-Kp2
    a32=Kp3;a52=-Kp3-Kp4;a62=Kp4
    a43=Kp3;a53=Kp4;a63=-Kp3-Kp4 
    #print(round(F0),round(F1,2),round(F2,2),round(F3,2),round(F4,2),round(F5,2),round(F6,2),round(F7,2))
    for xx1 in range(2):
        ta1=-1
        for xck in fen0:
            ta1+=1
            if xx1==0: 
                if ta1<len(dzk[0]):
                    cnnk=monum-num3-fen0[-ta1-1]
                if ta1>=len(dzk[0]):
                    break
            if xx1==1:
                if ta1<len(dzk[1]):
                    cnnk=xck-num4
                if ta1>=len(dzk[1]):
                    break
            for xx2 in range(2):
                ta2=-1
                for xckp in fenp0:
                    ta2+=1
                    if xx2==0:
                        if ta2<len(dzkp[0]):
                            cnnkp=xckp-nump4
                        if ta2>=len(dzkp[0]):
                            break
                    if xx2==1:
                        if ta2<len(dzkp[1]):
                            cnnkp=monump-nump3-fenp0[-ta2-1]
                        if ta2>=len(dzkp[1]):
                            break
                    for xx3 in range(2):
                        ta5=-1
                        for xckpx in fenpx0:    
                            ta5+=1
                            if xx3==0:
                                if ta5<len(dzkpx[0]):
                                    cnnkpx=xckpx-numpx4
                                if ta5>=len(dzkpx[0]):
                                    break
                            if xx3==1:
                                if ta5<len(dzkpx[1]):
                                    cnnkpx=monumpx-numpx3-fenpx0[-ta5-1]
                                if ta5>=len(dzkpx[1]):
                                    break
                            if xx1==0:
                                xa3=fen2
                            if xx1==1:
                                xa3=fen1
                            ta3=-1
                            for xc0 in xa3:
                                #print(xqq,xc0,len(dz0[xqq][1]),len(dz0[xqq][0]),'3')
                                ta3+=1  
                                if xx1==0:
                                    if ta3<len(dz0[0]):
                                        cnn0=num4-fen2[-ta3-1]
                                    if ta3>=len(dz0[0]):
                                        break
                                if xx1==1:
                                    if ta3<len(dz0[1]):
                                        cnn0=xc0-monum+num3
                                    if ta3>=len(dz0[1]):
                                        break
                                if xx2==0:
                                    xa4=fenp1
                                if xx2==1:
                                    xa4=fenp2
                                ta4=-1
                                for xc0p in xa4:
                                    ta4+=1
                                    if xx2==0:
                                        if ta4<len(dz0p[0]):
                                            cnn0p=xc0p-monump+nump3
                                        if ta4>=len(dz0p[0]):
                                            break
                                    if xx2==1:
                                        if ta4<len(dz0p[1]):
                                            cnn0p=nump4-fenp2[-ta4-1]    
                                    if xx3==0:
                                        xa5=fenpx1
                                    if xx3==1:
                                        xa5=fenpx2
                                    ta6=-1
                                    for xc0px in xa5:
                                        ta6+=1
                                        if xx3==0:
                                            if ta6<len(dz0px[0]):
                                                cnn0px=xc0px-monumpx+numpx3
                                            if ta6>=len(dz0px[0]):
                                                break
                                        if xx3==1:
                                            if ta6<len(dz0px[1]):
                                                cnn0px=numpx4-fenpx2[-ta6-1] 
                                            if ta6>len(dz0px[1]):
                                                break
                                        a04=-Kp1;a24=round(-(nump3+nump4+cnnkp-cnn0p)*Km,3);a34=round((nump3+nump4+cnnkp-cnn0p)*Km+Kp1+Kp3,3);a54=-Kp3
                                        a05=-Kp2;a15=round(((numpx3+numpx4+cnnkpx-cnn0px)+(num3+num4+cnnk-cnn0))*Km+Kp2,3);
                                        a25=round(-Km*(num3+num4+cnnk-cnn0),3);a45=round(-(numpx3+numpx4+cnnkpx-cnn0px)*Km,3)
                                        a16=round((num3+num4+cnnk-cnn0)*Km,3);a26=round(-Km*(nump3+nump4+cnnkp-cnn0p+num3+num4+cnnk-cnn0)-Kp2,3)
                                        a36=round(Km*(nump3+nump4+cnnkp-cnn0p),3);a76=Kp2
                                        a17=round(Km*(numpx3+numpx4+cnnkpx-cnn0px),3);a47=round(-Km*(numpx3+numpx4+cnnkpx-cnn0px)-Kp1-Kp3,3);a67=Kp3;a77=Kp1
                                        b4=-Km*(dzkp[xx2][ta2]-dz0p[xx2][ta4])-F1
                                        b5=-Km*(-dzkpx[xx3][ta5]+dz0px[xx3][ta6]+dzk[xx1][ta1]-dz0[xx1][ta3])-F2
                                        b6=-Km*(dzkp[xx2][ta2]-dz0p[xx2][ta4]+dzk[xx1][ta1]-dz0[xx1][ta3])-F3
                                        b7=Km*(dzkpx[xx3][ta5]-dz0px[xx3][ta6])-F4
                                        b=np.array([round(-F0,14),round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14)])
                                           
                                        try:
                                            A=np.array([[a00,a10,0,a30,0,0,0,0],[0,0,a21,0,a41,0,0,a71],[0,0,0,a32,0,a52,a62,0],
                                                        [0,0,0,0,a43,a53,a63,0],[a04,0,a24,a34,0,a54,0,0],[a05,a15,a25,0,a45,0,0,0],
                                                        [0,a16,a26,a36,0,0,0,a76],[0,a17,0,0,a47,0,a67,a77]])
                                             #b=np.array([-F0,-F7,-F5,-F6,b4,b5,b6,b7])

                                            mtui=np.linalg.solve(A,b)
                                        except:
                                            A=np.array([[0,a21,0,a41,0,0,a71],[0,0,a32,0,a52,a62,0],
                                                    [0,0,0,a43,a53,a63,0],[0,a24,a34,0,a54,0,0],[a15,a25,0,a45,0,0,0],
                                                    [a16,a26,a36,0,0,0,a76],[a17,0,0,a47,0,a67,a77]])
                                            b=np.array([round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14)])
                                          
                                            xtui=np.linalg.solve(A,b)
                                            for xn in range(7):
                                                mtui[xn+1]=xtui[xn]
                                            mtui[0]=0
                                        tuzx=(mtui[0]+mtui[7])/2
                                        for xn in range(8):
                                            mtui[xn]-=tuzx

                                        
                                        
                                        check0=-1#################
                                        checkk=-1################
                                        checkp0=-1#################
                                        checkpk=-1#################
                                        checkpx0=-1#################
                                        checkpxk=-1#################
                                        if xx1==0 and mtui[1]-mtui[2]>=0:#
                                            check0=0#################
                                            checkk=0#################
                                            for ton1 in range(monum):
                                                if jisuan[ton1]<-xdmotor and abs(jisuan[ton1]+mtui[1]-mtui[2])<=xdmotor:
                                                    check0+=1
                                                if abs(jisuan[ton1])<=xdmotor and jisuan[ton1]+mtui[1]-mtui[2]>xdmotor:
                                                    checkk+=1
                                        if xx1==1 and mtui[1]-mtui[2]<=0:
                                            check0=0#################
                                            checkk=0#################
                                            for ton1 in range(monum):
                                                if abs(jisuan[ton1])<=xdmotor and jisuan[ton1]+mtui[1]-mtui[2]<-xdmotor:
                                                    checkk+=1
                                                if jisuan[ton1]>xdmotor and abs(jisuan[ton1]+mtui[1]-mtui[2])<=xdmotor:
                                                    check0+=1
                                        if check0!=cnn0 or checkk!=cnnk:
                                            continue
                                        if xx2==1 and mtui[3]-mtui[2]>=0:
                                            checkp0=0#################
                                            checkpk=0#################
                                            for ton1 in range(monump):
                                                if jisuanp[ton1]<-xdmotor and abs(jisuanp[ton1]+mtui[3]-mtui[2])<=xdmotor:
                                                    checkp0+=1
                                                if abs(jisuanp[ton1])<=xdmotor and jisuanp[ton1]+mtui[3]-mtui[2]>xdmotor:
                                                    checkpk+=1
                                        if xx2==0 and mtui[3]-mtui[2]<=0:#
                                            checkp0=0#################
                                            checkpk=0#################
                                            for ton1 in range(monump):
                                                if abs(jisuanp[ton1])<=xdmotor and jisuanp[ton1]+mtui[3]-mtui[2]<-xdmotor:
                                                    checkpk+=1
                                                if jisuanp[ton1]>xdmotor and abs(jisuanp[ton1]+mtui[3]-mtui[2])<=xdmotor:
                                                    checkp0+=1
                                        if checkp0!=cnn0p or checkpk!=cnnkp:
                                            continue      
                                        if xx3==1 and mtui[4]-mtui[1]>=0:
                                            checkpx0=0################
                                            checkpxk=0################
                                            for ton1 in range(monumpx):
                                                if jisuanpx[ton1]<-xdmotor and abs(jisuanpx[ton1]+mtui[4]-mtui[1])<=xdmotor:
                                                    checkpx0+=1
                                                if abs(jisuanpx[ton1])<=xdmotor and jisuanpx[ton1]+mtui[4]-mtui[1]>xdmotor:
                                                    checkpxk+=1
                                        if xx3==0 and mtui[4]-mtui[1]<=0:#
                                            checkpx0=0#################
                                            checkpxk=0#################
                                            for ton1 in range(monumpx):
                                                if abs(jisuanpx[ton1])<=xdmotor and jisuanpx[ton1]+mtui[4]-mtui[1]<-xdmotor:
                                                    checkpxk+=1
                                                if jisuanpx[ton1]>xdmotor and abs(jisuanpx[ton1]+mtui[4]-mtui[1])<=xdmotor:
                                                    checkpx0+=1
                                        if checkpx0!=cnn0px or checkpxk!=cnnkpx:
                                            continue  
                                        #print(mtui,'qqq')
                                        return 1,mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]
    return 0,mtui[0],mtui[1],mtui[2],mtui[3],mtui[4],mtui[5],mtui[6],mtui[7]




def energy(i):#calculate the energy of the current spindle system and the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by two iMTs taking a forward step and a backward step 
    ypsl0=Kp3*(((Xkin-l1r)**2+(Xkin1-l4l)**2)/2)+Kp2*((Xpole-l2l)**2+(Xpo-l3r)**2)/2+Kp1*((Xpole-l1l)**2+(Xpo-l4r)**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
    for xn in range(nu):
        if kinb[xn]==1 and kina[xn]==1:
            dxc[xn]=(motor[xn][0]-motor[xn][1])
            if dxc[xn]>xdmotor:
                ypsl0+=Km*((dxc[xn]-xdmotor)**2)/2
            if dxc[xn]<-xdmotor:
                ypsl0+=Km*((dxc[xn]+xdmotor)**2)/2
        if kina[xn]==0 or kinb[xn]==0:
            dxc[xn]=0
    for xn in range(nu1):
        if kinbp[xn]==1 and kinap[xn]==1:
            if motorp[xn][0]-motorp[xn][1]>xdmotor:
                ypsl0+=Km*((motorp[xn][0]-motorp[xn][1]-xdmotor)**2)/2
            if motorp[xn][0]-motorp[xn][1]<-xdmotor:
                ypsl0+=Km*((motorp[xn][0]-motorp[xn][1]+xdmotor)**2)/2
    for xn in range(nu1x):
        if kinbpx[xn]==1 and kinapx[xn]==1:
            if motorpx[xn][0]-motorpx[xn][1]>xdmotor:
                ypsl0+=Km*((motorpx[xn][0]-motorpx[xn][1]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]<-xdmotor:
                ypsl0+=Km*((motorpx[xn][0]-motorpx[xn][1]+xdmotor)**2)/2
    ypsl11=Kp3*(((Xkin+zj5[i]-l1r-zj3[i])**2+(Xkin1+zj6[i]-l4l-zj4[i])**2)/2)+Kp2*((Xpole+zj0[i]-l2l-zj1[i])**2+(Xpo+zj7[i]-l3r-zj2[i])**2)/2+Kp1*((Xpole+zj0[i]-l1l-zj3[i])**2+(Xpo+zj7[i]-l4r-zj4[i])**2)/2+Kp4*((Xkin1+zj6[i]-Xkin-zj5[i]-kdis)**2)/2
    ypsl22=Kp3*(((Xkin+js5[i]-l1r-js3[i])**2+(Xkin1+js6[i]-l4l-js4[i])**2)/2)+Kp2*((Xpole+js0[i]-l2l-js1[i])**2+(Xpo+js7[i]-l3r-js2[i])**2)/2+Kp1*((Xpole+js0[i]-l1l-js3[i])**2+(Xpo+js7[i]-l4r-js4[i])**2)/2+Kp4*((Xkin1+js6[i]-Xkin-js5[i]-kdis)**2)/2
    
    for xn in range(nu):
        if kina[xn]==1 and kinb[xn]==1:
            if xn==i:
                if dx1+d+zj1[i]-zj2[i]>xdmotor:
                    ypsl11+=Km*((dx1+d+zj1[i]-zj2[i]-xdmotor)**2)/2
                if dx1+d+zj1[i]-zj2[i]<-xdmotor:
                    ypsl11+=Km*((dx1+d+zj1[i]-zj2[i]+xdmotor)**2)/2
                if dx1-d+js1[i]-js2[i]>xdmotor:
                    ypsl22+=Km*((dx1-d+js1[i]-js2[i]-xdmotor)**2)/2
                if dx1-d+js1[i]-js2[i]<-xdmotor:
                    ypsl22+=Km*((dx1-d+js1[i]-js2[i]+xdmotor)**2)/2
            if xn!=i:
                if dxc[xn]+zj1[i]-zj2[i]>xdmotor:
                    ypsl11+=Km*((dxc[xn]+zj1[i]-zj2[i]-xdmotor)**2)/2
                if dxc[xn]+zj1[i]-zj2[i]<-xdmotor:
                    ypsl11+=Km*((dxc[xn]+zj1[i]-zj2[i]+xdmotor)**2)/2
                if dxc[xn]+js1[i]-js2[i]>xdmotor:
                    ypsl22+=Km*((dxc[xn]+js1[i]-js2[i]-xdmotor)**2)/2
                if dxc[xn]+js1[i]-js2[i]<-xdmotor:
                    ypsl22+=Km*((dxc[xn]+js1[i]-js2[i]+xdmotor)**2)/2
    for xn in range(nu1):
        if kinap[xn]+kinbp[xn]==2:
            if motorp[xn][0]-motorp[xn][1]+zj3[i]-zj2[i]>xdmotor:
                ypsl11+=Km*((motorp[xn][0]-motorp[xn][1]+zj3[i]-zj2[i]-xdmotor)**2)/2
            if motorp[xn][0]-motorp[xn][1]+zj3[i]-zj2[i]<-xdmotor:
                ypsl11+=Km*((motorp[xn][0]-motorp[xn][1]+zj3[i]-zj2[i]+xdmotor)**2)/2
            if motorp[xn][0]-motorp[xn][1]+js3[i]-js2[i]>xdmotor:
                ypsl22+=Km*((motorp[xn][0]-motorp[xn][1]+js3[i]-js2[i]-xdmotor)**2)/2
            if motorp[xn][0]-motorp[xn][1]+js3[i]-js2[i]<-xdmotor:
                ypsl22+=Km*((motorp[xn][0]-motorp[xn][1]+js3[i]-js2[i]+xdmotor)**2)/2
    for xn in range(nu1x):
        if kinapx[xn]+kinbpx[xn]==2:
            if motorpx[xn][0]-motorpx[xn][1]+zj4[i]-zj1[i]>xdmotor:
                ypsl11+=Km*((motorpx[xn][0]-motorpx[xn][1]+zj4[i]-zj1[i]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+zj4[i]-zj1[i]<-xdmotor:
                ypsl11+=Km*((motorpx[xn][0]-motorpx[xn][1]+zj4[i]-zj1[i]+xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+js4[i]-js1[i]>xdmotor:
                ypsl22+=Km*((motorpx[xn][0]-motorpx[xn][1]+js4[i]-js1[i]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+js4[i]-js1[i]<-xdmotor:
                ypsl22+=Km*((motorpx[xn][0]-motorpx[xn][1]+js4[i]-js1[i]+xdmotor)**2)/2
    return ypsl0,ypsl11,ypsl22


# In[8]:


def energyp(i):#calculate the energy of the current spindle system and the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by the iMT connected to the right spindle pole and the kMT connected to the left spindle pole taking a forward step and a backward step 
    ypsl0p=Kp3*(((Xkin-l1r)**2+(Xkin1-l4l)**2)/2)+Kp2*((Xpole-l2l)**2+(Xpo-l3r)**2)/2+Kp1*((Xpole-l1l)**2+(Xpo-l4r)**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
    for xn in range(nu):
        if kinb[xn]==1 and kina[xn]==1:
            dxc[xn]=(motor[xn][0]-motor[xn][1])
            if dxc[xn]>xdmotor:
                ypsl0p+=Km*((dxc[xn]-xdmotor)**2)/2
            if dxc[xn]<-xdmotor:
                ypsl0p+=Km*((dxc[xn]+xdmotor)**2)/2
        if kina[xn]==0 or kinb[xn]==0:
            dxc[xn]=0
    for xn in range(nu1):
        if kinbp[xn]==1 and kinap[xn]==1:
            if motorp[xn][0]-motorp[xn][1]>xdmotor:
                ypsl0p+=Km*((motorp[xn][0]-motorp[xn][1]-xdmotor)**2)/2
            if motorp[xn][0]-motorp[xn][1]<-xdmotor:
                ypsl0p+=Km*((motorp[xn][0]-motorp[xn][1]+xdmotor)**2)/2
    for xn in range(nu1x):
        if kinbpx[xn]==1 and kinapx[xn]==1:
            if motorpx[xn][0]-motorpx[xn][1]>xdmotor:
                ypsl0p+=Km*((motorpx[xn][0]-motorpx[xn][1]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]<-xdmotor:
                ypsl0p+=Km*((motorpx[xn][0]-motorpx[xn][1]+xdmotor)**2)/2
    #ypsl1p1=Kp3*(((Xkin+zj5p[i]-l1r-zj3p[i])**2+(Xkin1+zj6p[i]-l4l-zj4p[i])**2)/2)+Kp2*((Xpole+zj0p[i]-l2l-zj1p[i])**2+(Xpo+zj7p[i]-l3r-zj2p[i])**2)/2+Kp1*((Xpole+zj0p[i]-l1l-zj3p[i])**2+(Xpo+zj7p[i]-l4r-zj4p[i])**2)/2+Kp4*((Xkin1+zj6p[i]-Xkin-zj5p[i]-kdis)**2)/2
    #ypsl2p2=Kp3*(((Xkin+js5p[i]-l1r-js3p[i])**2+(Xkin1+js6p[i]-l4l-js4p[i])**2)/2)+Kp2*((Xpole+js0p[i]-l2l-js1p[i])**2+(Xpo+js7p[i]-l3r-js2p[i])**2)/2+Kp1*((Xpole+js0p[i]-l1l-js3p[i])**2+(Xpo+js7p[i]-l4r-js4p[i])**2)/2+Kp4*((Xkin1+js6p[i]-Xkin-js5p[i]-kdis)**2)/2
    
    ypsl1p1=Kp3*(((Xkin+zj5p[i]-l1r-zj3p[i])**2+(Xkin1+zj6p[i]-l4l-zj4p[i])**2)/2)+Kp2*((Xpole+zj0p[i]-l2l-zj1p[i])**2+(Xpo+zj7p[i]-l3r-zj2p[i])**2)/2+Kp1*((Xpole+zj0p[i]-l1l-zj3p[i])**2+(Xpo+zj7p[i]-l4r-zj4p[i])**2)/2+Kp4*((Xkin1+zj6p[i]-Xkin-zj5p[i]-kdis)**2)/2
    ypsl2p2=Kp3*(((Xkin+js5p[i]-l1r-js3p[i])**2+(Xkin1+js6p[i]-l4l-js4p[i])**2)/2)+Kp2*((Xpole+js0p[i]-l2l-js1p[i])**2+(Xpo+js7p[i]-l3r-js2p[i])**2)/2+Kp1*((Xpole+js0p[i]-l1l-js3p[i])**2+(Xpo+js7p[i]-l4r-js4p[i])**2)/2+Kp4*((Xkin1+js6p[i]-Xkin-js5p[i]-kdis)**2)/2
    #print(ypsl1p1,ypsl2p2)
    for xn in range(nu):
        if kina[xn]==1 and kinb[xn]==1:
            if dxc[xn]+zj1p[i]-zj2p[i]>xdmotor:
                ypsl1p1+=Km*((dxc[xn]+zj1p[i]-zj2p[i]-xdmotor)**2)/2
            if dxc[xn]+zj1p[i]-zj2p[i]<-xdmotor:
                ypsl1p1+=Km*((dxc[xn]+zj1p[i]-zj2p[i]+xdmotor)**2)/2
            if dxc[xn]+js1p[i]-js2p[i]>xdmotor:
                ypsl2p2+=Km*((dxc[xn]+js1p[i]-js2p[i]-xdmotor)**2)/2
            if dxc[xn]+js1p[i]-js2p[i]<-xdmotor:
                ypsl2p2+=Km*((dxc[xn]+js1p[i]-js2p[i]+xdmotor)**2)/2
    for xn in range(nu1):
        if kinap[xn]+kinbp[xn]==2:
            if xn==i:
                if dx1+d+zj3p[i]-zj2p[i]>xdmotor:
                    ypsl1p1+=Km*((dx1+d+zj3p[i]-zj2p[i]-xdmotor)**2)/2
                if dx1+d+zj3p[i]-zj2p[i]<-xdmotor:
                    ypsl1p1+=Km*((dx1+d+zj3p[i]-zj2p[i]+xdmotor)**2)/2
                if dx1-d+js3p[i]-js2p[i]>xdmotor:
                    ypsl2p2+=Km*((dx1-d+js3p[i]-js2p[i]-xdmotor)**2)/2
                if dx1-d+js3p[i]-js2p[i]<-xdmotor:
                    ypsl2p2+=Km*((dx1-d+js3p[i]-js2p[i]+xdmotor)**2)/2
            if xn!=i:    
                if motorp[xn][0]-motorp[xn][1]+zj3p[i]-zj2p[i]>xdmotor:
                    ypsl1p1+=Km*((motorp[xn][0]-motorp[xn][1]+zj3p[i]-zj2p[i]-xdmotor)**2)/2
                if motorp[xn][0]-motorp[xn][1]+zj3p[i]-zj2p[i]<-xdmotor:
                    ypsl1p1+=Km*((motorp[xn][0]-motorp[xn][1]+zj3p[i]-zj2p[i]+xdmotor)**2)/2
                if motorp[xn][0]-motorp[xn][1]+js3p[i]-js2p[i]>xdmotor:
                    ypsl2p2+=Km*((motorp[xn][0]-motorp[xn][1]+js3p[i]-js2p[i]-xdmotor)**2)/2
                if motorp[xn][0]-motorp[xn][1]+js3p[i]-js2p[i]<-xdmotor:
                    ypsl2p2+=Km*((motorp[xn][0]-motorp[xn][1]+js3p[i]-js2p[i]+xdmotor)**2)/2
    for xn in range(nu1x):
        if kinapx[xn]+kinbpx[xn]==2:
            if motorpx[xn][0]-motorpx[xn][1]+zj4p[i]-zj1p[i]>xdmotor:
                ypsl1p1+=Km*((motorpx[xn][0]-motorpx[xn][1]+zj4p[i]-zj1p[i]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+zj4p[i]-zj1p[i]<-xdmotor:
                ypsl1p1+=Km*((motorpx[xn][0]-motorpx[xn][1]+zj4p[i]-zj1p[i]+xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+js4p[i]-js1p[i]>xdmotor:
                ypsl2p2+=Km*((motorpx[xn][0]-motorpx[xn][1]+js4p[i]-js1p[i]-xdmotor)**2)/2
            if motorpx[xn][0]-motorpx[xn][1]+js4p[i]-js1p[i]<-xdmotor:
                ypsl2p2+=Km*((motorpx[xn][0]-motorpx[xn][1]+js4p[i]-js1p[i]+xdmotor)**2)/2
    return ypsl0p,ypsl1p1,ypsl2p2


# In[9]:


def energypx(i):#calculate the energy of the current spindle system and the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by the iMT connected to the left spindle pole and the kMT connected to the right spindle pole taking a forward step and a backward step 
    ypsl0p=Kp3*(((Xkin-l1r)**2+(Xkin1-l4l)**2)/2)+Kp2*((Xpole-l2l)**2+(Xpo-l3r)**2)/2+Kp1*((Xpole-l1l)**2+(Xpo-l4r)**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
    for xn in range(nu):
        if kinb[xn]==1 and kina[xn]==1:
            dxc[xn]=(motor[xn][0]-motor[xn][1])
            if dxc[xn]>xdmotor:
                ypsl0p+=Km*((dxc[xn]-xdmotor)**2)/2
            if dxc[xn]<-xdmotor:
                ypsl0p+=Km*((dxc[xn]+xdmotor)**2)/2
        if kina[xn]==0 or kinb[xn]==0:
            dxc[xn]=0
    for xn in range(nu1):
        if kinbp[xn]==1 and kinap[xn]==1:
            dxcp[xn]=motorp[xn][0]-motorp[xn][1]
            if dxcp[xn]>xdmotor:
                ypsl0p+=Km*((dxcp[xn]-xdmotor)**2)/2
            if dxcp[xn]<-xdmotor:
                ypsl0p+=Km*((dxcp[xn]+xdmotor)**2)/2
        if kinbp[xn]==0 or kinap[xn]==0:
            dxcp[xn]=0
    for xn in range(nu1x):
        if kinbpx[xn]==1 and kinapx[xn]==1:
            dxcpx[xn]=motorpx[xn][0]-motorpx[xn][1]
            if dxcpx[xn]>xdmotor:
                ypsl0p+=Km*((dxcpx[xn]-xdmotor)**2)/2
            if dxcpx[xn]<-xdmotor:
                ypsl0p+=Km*((dxcpx[xn]+xdmotor)**2)/2
        if kinbpx[xn]==0 or kinapx[xn]==0:
            dxcpx[xn]=0
    ypsl1p1=Kp3*(((Xkin+zj5px[i]-l1r-zj3px[i])**2+(Xkin1+zj6px[i]-l4l-zj4px[i])**2)/2)+Kp2*((Xpole+zj0px[i]-l2l-zj1px[i])**2+(Xpo+zj7px[i]-l3r-zj2px[i])**2)/2+Kp1*((Xpole+zj0px[i]-l1l-zj3px[i])**2+(Xpo+zj7px[i]-l4r-zj4px[i])**2)/2+Kp4*((Xkin1+zj6px[i]-Xkin-zj5px[i]-kdis)**2)/2
    ypsl2p2=Kp3*(((Xkin+js5px[i]-l1r-js3px[i])**2+(Xkin1+js6px[i]-l4l-js4px[i])**2)/2)+Kp2*((Xpole+js0px[i]-l2l-js1px[i])**2+(Xpo+js7px[i]-l3r-js2px[i])**2)/2+Kp1*((Xpole+js0px[i]-l1l-js3px[i])**2+(Xpo+js7px[i]-l4r-js4px[i])**2)/2+Kp4*((Xkin1+js6px[i]-Xkin-js5px[i]-kdis)**2)/2
    for xn in range(nu):
        if kina[xn]==1 and kinb[xn]==1:
            if dxc[xn]+zj1px[i]-zj2px[i]>xdmotor:
                ypsl1p1+=Km*((dxc[xn]+zj1px[i]-zj2px[i]-xdmotor)**2)/2
            if dxc[xn]+zj1px[i]-zj2px[i]<-xdmotor:
                ypsl1p1+=Km*((dxc[xn]+zj1px[i]-zj2px[i]+xdmotor)**2)/2
            if dxc[xn]+js1px[i]-js2px[i]>xdmotor:
                ypsl2p2+=Km*((dxc[xn]+js1px[i]-js2px[i]-xdmotor)**2)/2
            if dxc[xn]+js1px[i]-js2px[i]<-xdmotor:
                ypsl2p2+=Km*((dxc[xn]+js1px[i]-js2px[i]+xdmotor)**2)/2
    for xn in range(nu1):
        if kinap[xn]+kinbp[xn]==2:
            if dxcp[xn]+zj3px[i]-zj2px[i]>xdmotor:
                ypsl1p1+=Km*((dxcp[xn]+zj3px[i]-zj2px[i]-xdmotor)**2)/2
            if dxcp[xn]+zj3px[i]-zj2px[i]<-xdmotor:
                ypsl1p1+=Km*((dxcp[xn]+zj3px[i]-zj2px[i]+xdmotor)**2)/2
            if dxcp[xn]+js3px[i]-js2px[i]>xdmotor:
                ypsl2p2+=Km*((dxcp[xn]+js3px[i]-js2px[i]-xdmotor)**2)/2
            if dxcp[xn]+js3px[i]-js2px[i]<-xdmotor:
                ypsl2p2+=Km*((dxcp[xn]+js3px[i]-js2px[i]+xdmotor)**2)/2
    for xn in range(nu1x):
        
        if kinapx[xn]+kinbpx[xn]==2:
            if xn==i:
                if dxcpx[xn]+d+zj4px[i]-zj1px[i]>xdmotor:
                    ypsl1p1+=Km*((dxcpx[xn]+d+zj4px[i]-zj1px[i]-xdmotor)**2)/2
                if dxcpx[xn]+d+zj4px[i]-zj1px[i]<-xdmotor:
                    ypsl1p1+=Km*((dxcpx[xn]+d+zj4px[i]-zj1px[i]+xdmotor)**2)/2
                if dxcpx[xn]-d+js4px[i]-js1px[i]>xdmotor:
                    ypsl2p2+=Km*((dxcpx[xn]-d+js4px[i]-js1px[i]-xdmotor)**2)/2
                if dxcpx[xn]-d+js4px[i]-js1px[i]<-xdmotor:
                    ypsl2p2+=Km*((dxcpx[xn]-d+js4px[i]-js1px[i]+xdmotor)**2)/2
            if xn!=i:
                if dxcpx[xn]+zj4px[i]-zj1px[i]>xdmotor:
                    ypsl1p1+=Km*((dxcpx[xn]+zj4px[i]-zj1px[i]-xdmotor)**2)/2
                if dxcpx[xn]+zj4px[i]-zj1px[i]<-xdmotor:
                    ypsl1p1+=Km*((dxcpx[xn]+zj4px[i]-zj1px[i]+xdmotor)**2)/2
                if dxcpx[xn]+js4px[i]-js1px[i]>xdmotor:
                    ypsl2p2+=Km*((dxcpx[xn]+js4px[i]-js1px[i]-xdmotor)**2)/2
                if dxcpx[xn]+js4px[i]-js1px[i]<-xdmotor:
                    ypsl2p2+=Km*((dxcpx[xn]+js4px[i]-js1px[i]+xdmotor)**2)/2
    return ypsl0p,ypsl1p1,ypsl2p2


# In[10]:


def checkF():#Calculate the force of the two kinetochores, two spindle poles and four MTs
    F0=Kp1*(l1l-Xpole)+Kp2*(l2l-Xpole)#the force of the left spindle pole
    F7=Kp1*(l4r-Xpo)+Kp2*(l3r-Xpo)#the force of the right spindle pole
    F5=Kp3*(l1r-Xkin)+Kp4*(Xkin1-Xkin-kdis)#the force of the left kinetochore
    F6=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l-Xkin1)#the force of the right kinetochore
    F3=Kp2*(Xpo-l3r)#the force of the iMT connected to the right spindle pole
    F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)#the force of the kMT connected to the left spindle pole
    F2=Kp2*(l2l-Xpole)#the force of the iMT connected to the left spindle pole
    F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)#the force of the kMT connected to the right spindle pole
    for ton1 in range(nu):
        if kina[ton1]==1 and kinb[ton1]==1:
            dxc[ton1]=motor[ton1][0]-motor[ton1][1]
            if dxc[ton1]>xdmotor:
                F3+=(dxc[ton1]-xdmotor)*Km
                F2+=(dxc[ton1]-xdmotor)*Km
            if dxc[ton1]<-xdmotor:
                F3+=(dxc[ton1]+xdmotor)*Km
                F2+=(dxc[ton1]+xdmotor)*Km
        if kina[ton1]==0 or kinb[ton1]==0:
            dxc[ton1]=0
    for ton1 in range(nu1):
        if kinap[ton1]==1 and kinbp[ton1]==1:
            dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
            if dxcp[ton1]>xdmotor:
                F3+=(dxcp[ton1]-xdmotor)*Km
                F1+=(dxcp[ton1]-xdmotor)*Km
            if dxcp[ton1]<-xdmotor:
                F3+=(dxcp[ton1]+xdmotor)*Km
                F1+=(dxcp[ton1]+xdmotor)*Km
        if kinap[ton1]==0 or kinbp[ton1]==0:
            dxcp[ton1]=0
    for ton1 in range(nu1x):
        if kinapx[ton1]==1 and kinbpx[ton1]==1:
            dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
            if dxcpx[ton1]>xdmotor:
                F4-=(dxcpx[ton1]-xdmotor)*Km
                F2-=(dxcpx[ton1]-xdmotor)*Km
            if dxcpx[ton1]<-xdmotor:
                F4-=(dxcpx[ton1]+xdmotor)*Km
                F2-=(dxcpx[ton1]+xdmotor)*Km
        if kinapx[ton1]==0 or kinbpx[ton1]==0:
            dxcpx[ton1]=0
    if abs(F1)+abs(F2)+abs(F3)+abs(F4)+abs(F0)+abs(F7)+abs(F5)+abs(F6)>0.1:
        de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
        return 1,de0,de1,de2,de3,de4,de5,de6,de7
    return 0,0,0,0,0,0,0,0,0


# main

##For storing intermediate data
MTk1=[]
MTi1=[]
MTzn1=[]
MTzn2=[]
MTzn3=[]
MTzn4=[]
MTzn5=[]
MTzn6=[]
MTzn7=[]
MTzn8=[]
MTzn9=[]
MTzn10=[]
MTzn11=[]
MTzn12=[]
MTzn13=[]
MTzn14=[]
MTzn15=[]
time1=[]
MTX1=[]
MTX2=[]	


h = 0.001
Km = 0.55# pN/nm
d=8.2
dcc=8.2
Number=5
vp1=4
vpol=vp1/d
xdmotor=d+0.000001
nant=3#[Eg5]
npar=1
ka1=0.0004*nant
ka2=0.0004*nant
Kp2=0.1
Kp1=0.1
Kp3=0.1
Kp4=10
lover=5002
kdis=1000.4
llap=4001.6
l1l=-lover#the kMT connected to the left spindle pole
l1r=-kdis/2
l2l=-lover#the iMT connected to the left spindle pole
l2r=llap
l3l=-llap#the iMT connected to the right spindle pole
l3r=-l2l
l4l=kdis/2#the kMT connected to the right spindle pole
l4r=l3r
Xpole=l1l
Xpo=l4r
Xkin=l1r
Xkin1=l4l
pbei=4#B
vp0=vp1/pbei
Fp0=3.5
#parameters of the MCAK
kon=0.0001
koff=0.001
kdif=82
vi=5
vk=4
stai=0.5
stad=0.5
##################
Ximt2=0
Ximt3=0
Xkmt1=0
Xkmt4=0
desi2=0
desi3=0
desk1=0
desk4=0
pos1=0
pos4=0
desia2=0
desia3=0
deska1=0
deska4=0
posa1=0
posa4=0
ximta2=0
ximta3=0
xkmta1=0
xkmta4=0
ads=0
fenbu=[]
adeps1=0
adeps2=0
adeps3=0
adeps4=0
panx=0
site=round(lover/d)-1
kina =[]
kinb =[]
kinap =[]
kinbp =[]
kinapx=[]
kinbpx=[]
motor =[]#kinesin-5 in the overlap region formed by two iMTs
motorp=[]#kinesin-5 in the overlap region formed by the iMT connected to the right spindle pole and the kMT connected to the left spindle pole 
motorpx=[]#the kinesin-5 motors in the overlap region formed by the iMT connected to the left spindle pole and the kMT connected to the right spindle pole 
motormc1=[]
motormc2=[]
motormc3=[]
motormc4=[]
dxc =[]
dxcp =[]
dxcpx =[]
pan=[]
ypsl1=[]
ypsl2=[]
panp=[]
ypsl1p=[]
ypsl2p=[]
panpx=[]
ypsl1px=[]
ypsl2px=[]
zj0=[]
js0=[]
zj1=[]
js1=[]
zj2=[]
js2=[]
zj3=[]
js3=[]
zj4=[]
js4=[]
zj5=[]
js5=[]
zj6=[]
js6=[]
zj7=[]
js7=[]
zj0p=[]
js0p=[]
zj1p=[]
js1p=[]
zj2p=[]
js2p=[]
zj3p=[]
js3p=[]
zj4p=[]
js4p=[]
zj5p=[]
js5p=[]
zj6p=[]
js6p=[]
zj7p=[]
js7p=[]
zj0px=[]
js0px=[]
zj1px=[]
js1px=[]
zj2px=[]
js2px=[]
zj3px=[]
js3px=[]
zj4px=[]
js4px=[]
zj5px=[]
js5px=[]
zj6px=[]
js6px=[]
zj7px=[]
js7px=[]
xx=1/8
for xn in range(4*7):
    xx+=1/4/8
    cx0=l3r-round((lover-llap)/d*random.uniform(0,xx))*d
    motormc3.append(cx0)
xx=1/8
for xn in range(4*7):
    xx+=1/4/8
    cx0=l1l+round((lover-llap)/d*random.uniform(0,xx))*d
    motormc1.append(cx0)
xx=1/8
for xn in range(4*7):
    xx+=1/4/8
    cx0=l4r-round((lover-llap)/d*random.uniform(0,xx))*d
    motormc4.append(cx0)
xx=1/8
for xn in range(4*7):
    xx+=1/4/8
    cx0=l2l+round((lover-llap)/d*random.uniform(0,xx))*d
    motormc2.append(cx0)
for xn in range(round(Number*random.uniform(nant*2,nant*2+1))):
    cx0=l3l+round(lover/d*random.uniform(0,1))*d
    motor.append([cx0,cx0])
    kina.append(round(random.uniform(0,1)))
    kinb.append(round(random.uniform(0,1)))
    ypsl1.append(0)
    dxc.append(0)
    pan.append(0)
    ypsl2.append(0)
    zj0.append(0)
    js0.append(0)
    zj1.append(0)
    js1.append(0)
    zj2.append(0)
    js2.append(0)
    zj3.append(0)
    js3.append(0)
    zj4.append(0)
    js4.append(0)
    zj5.append(0)
    js5.append(0)
    zj6.append(0)
    js6.append(0)
    zj7.append(0)
    js7.append(0)

for xn in range(round(Number*random.uniform(nant,nant+1))):
    cx0=l3l+round(lover/2/d*random.uniform(0,1))*d
    motorp.append([cx0,cx0])
    kinap.append(round(random.uniform(0,1)))
    kinbp.append(round(random.uniform(0,1)))
    dxcp.append(0)
    panp.append(0)
    ypsl1p.append(0)
    ypsl2p.append(0)
    zj0p.append(0)
    js0p.append(0)
    zj1p.append(0)
    js1p.append(0)
    zj2p.append(0)
    js2p.append(0)
    zj3p.append(0)
    js3p.append(0)
    zj4p.append(0)
    js4p.append(0)
    zj5p.append(0)
    js5p.append(0)
    zj6p.append(0)
    js6p.append(0)
    zj7p.append(0)
    js7p.append(0)

for xn in range(round(Number*random.uniform(nant,nant+1))):
    cx0=l4l+round(lover/2/d*random.uniform(0,1))*d
    motorpx.append([cx0,cx0])
    kinapx.append(round(random.uniform(0,1)))
    kinbpx.append(round(random.uniform(0,1)))
    dxcpx.append(0)
    panpx.append(0)
    ypsl1px.append(0)
    ypsl2px.append(0)
    zj0px.append(0)
    js0px.append(0)
    zj1px.append(0)
    js1px.append(0)
    zj2px.append(0)
    js2px.append(0)
    zj3px.append(0)
    js3px.append(0)
    zj4px.append(0)
    js4px.append(0)
    zj5px.append(0)
    js5px.append(0)
    zj6px.append(0)
    js6px.append(0)
    zj7px.append(0)
    js7px.append(0)
nu=len(kina)
nu1=len(kinap)
nu1x=len(kinapx)
pan1=1
monum=nu
monump=nu1
monumpx=nu1x


#loop


for a in range(0,10000000000):
#################################################################
    nu=len(kina)
    nu1=len(kinap)
    nu1x=len(kinapx)

    lspa1=min(l3l-l1l,Xkin-l1l)##1MT
    lspa2=l3l-l2l##2MT
    lspa3=l3r-l2r#3MT
    lspa4=min(l4r-l2r,l4r-Xkin1)
    numc1=len(motormc1)
    numc2=len(motormc2)
    numc3=len(motormc3)
    numc4=len(motormc4)
    nn1=0
    nn2=0
    nn3=0
    for xn in range(numc1):
        if motormc1[xn]<l1l+1000:
            nn1+=1
        if motormc1[xn]>=l1l+1000 and motormc1[xn]<l1l+2000:
            nn2+=1
        if motormc1[xn]>=l1l+2000 and motormc1[xn]<l1l+3000:
            nn3+=1
    if lspa1>1000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn1*d)*h:
            cx0=l1l+round(1000/d*random.uniform(0,1))*d
            it=0
            while it <numc1:
                if abs(cx0-motormc1[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc1.append(cx0) 
    if lspa1>2000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn2*d)*h:
            cx0=l1l+round(2000/d*random.uniform(0.5,1))*d
            it=0
            while it <numc1:
                if abs(cx0-motormc1[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc1.append(cx0) 
    if lspa1>3000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn3*d)*h:
            cx0=l1l+round(3000/d*random.uniform(0.67,1))*d
            it=0
            while it <numc1:
                if abs(cx0-motormc1[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc1.append(cx0)
    if lspa1<=1000:
        nc1=numc1
        ld=lspa1
    if lspa1<=2000 and lspa1>1000:
        nc1=numc1-nn1
        ld=lspa1-1000
    if lspa1<=3000 and lspa1>2000:
        nc1=numc1-nn1-nn2
        ld=lspa1-2000
    if lspa1>3000:
        nc1=numc1-nn1-nn2-nn3
        ld=lspa1-3000
    ran5 = random.uniform(0, 1)
    if ran5<kon*(ld-nc1*d)*h:
        cx0=l1l+round(lspa1/d*random.uniform(1-ld/lspa1,1))*d
        it=0
        while it <numc1:
            if abs(cx0-motormc1[it])<0:
                cx0+=d
                it=-1
            it+=1
        motormc1.append(cx0)
    nn1=0
    nn2=0
    nn3=0
    for xn in range(numc2):
        if motormc2[xn]<l2l+1000:
            nn1+=1
        if motormc2[xn]>=l2l+1000 and motormc2[xn]<l2l+2000:
            nn2+=1
        if motormc2[xn]>=l2l+2000 and motormc2[xn]<l2l+3000:
            nn3+=1
    if lspa2>1000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn1*d)*h:
            cx0=l2l+round(1000/d*random.uniform(0,1))*d
            it=0
            while it <numc2:
                if abs(cx0-motormc2[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc2.append(cx0) 
    if lspa2>2000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn2*d)*h:
            cx0=l2l+round(2000/d*random.uniform(0.5,1))*d
            it=0
            while it <numc2:
                if abs(cx0-motormc2[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc2.append(cx0) 
    if lspa2>3000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn3*d)*h:
            cx0=l2l+round(3000/d*random.uniform(0.67,1))*d
            it=0
            while it <numc2:
                if abs(cx0-motormc2[it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc2.append(cx0)
    if lspa2<=1000:
        nc1=numc2
        ld=lspa2
    if lspa2<=2000 and lspa2>1000:
        nc1=numc2-nn1
        ld=lspa2-1000
    if lspa2<=3000 and lspa2>2000:
        nc1=numc2-nn1-nn2
        ld=lspa2-2000
    if lspa2>3000:
        nc1=numc2-nn1-nn2-nn3
        ld=lspa2-3000
    ran5 = random.uniform(0, 1)
    if ran5<kon*(ld-nc1*d)*h:
        cx0=l2l+round(lspa2/d*random.uniform(1-ld/lspa2,1))*d
        it=0
        while it <numc2:
            if abs(cx0-motormc2[it])<0:
                cx0+=d
                it=-1
            it+=1
        motormc2.append(cx0)
    nn1=0
    nn2=0
    nn3=0
    for xn in range(numc3):
        if motormc3[xn]>l3r-1000:
            nn1+=1
        if motormc3[xn]<=l3r-1000 and motormc3[xn]>l3r-2000:
            nn2+=1
        if motormc3[xn]<=l3r-2000 and motormc3[xn]>l3r-3000:
            nn3+=1
    if lspa3>1000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn1*d)*h:
            cx0=l3r-round(1000/d*random.uniform(0,1))*d
            it=0
            while it <numc3:
                if abs(cx0-motormc3[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc3.append(cx0) 
    if lspa3>2000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn2*d)*h:
            cx0=l3r-round(2000/d*random.uniform(0.5,1))*d
            it=0
            while it <numc3:
                if abs(cx0-motormc3[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc3.append(cx0) 
    if lspa3>3000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn3*d)*h:
            cx0=l3r-round(3000/d*random.uniform(0.67,1))*d
            it=0
            while it <numc3:
                if abs(cx0-motormc3[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc3.append(cx0)
    if lspa3<=1000:
        nc1=numc3
        ld=lspa3
    if lspa3<=2000 and lspa3>1000:
        nc1=numc3-nn1
        ld=lspa3-1000
    if lspa3<=3000 and lspa3>2000:
        nc1=numc3-nn1-nn2
        ld=lspa3-2000
    if lspa3>3000:
        nc1=numc3-nn1-nn2-nn3
        ld=lspa3-3000
    ran4 = random.uniform(0, 1)
    if ran4<kon*(ld-nc1*d)*h:
        cx0=l3r-round(lspa3/d*random.uniform(1-ld/lspa3,1))*d
        it=0
        while it <numc3:
            if abs(cx0-motormc3[it])<0:
                cx0-=d
                it=-1
            it+=1
        motormc3.append(cx0)
    nn1=0
    nn2=0
    nn3=0
    for xn in range(numc4):
        if motormc4[xn]>l4r-1000:
            nn1+=1
        if motormc4[xn]<=l4r-1000 and motormc4[xn]>l4r-2000:
            nn2+=1
        if motormc4[xn]<=l4r-2000 and motormc4[xn]>l4r-3000:
            nn3+=1
    if lspa4>1000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn1*d)*h:
            cx0=l4r-round(1000/d*random.uniform(0,1))*d
            it=0
            while it <numc4:
                if abs(cx0-motormc4[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc4.append(cx0) 
    if lspa4>2000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn2*d)*h:
            cx0=l4r-round(2000/d*random.uniform(0.5,1))*d
            it=0
            while it <numc4:
                if abs(cx0-motormc4[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc4.append(cx0) 
    if lspa4>3000:
        ran5 = random.uniform(0, 1)
        if ran5<kon*(1000-nn3*d)*h:
            cx0=l4r-round(3000/d*random.uniform(0.67,1))*d
            it=0
            while it <numc4:
                if abs(cx0-motormc4[it])<0:
                    cx0-=d
                    it=-1
                it+=1
            motormc4.append(cx0)
    if lspa4<=1000:
        nc1=numc4
        ld=lspa4
    if lspa4<=2000 and lspa4>1000:
        nc1=numc4-nn1
        ld=lspa4-1000
    if lspa4<=3000 and lspa4>2000:
        nc1=numc4-nn1-nn2
        ld=lspa4-2000
    if lspa4>3000:
        nc1=numc4-nn1-nn2-nn3
        ld=lspa4-3000
    ran4 = random.uniform(0, 1)
    if ran4<kon*(ld-nc1*d)*h:
        cx0=l4r-round(lspa4/d*random.uniform(1-ld/lspa4,1))*d
        it=0
        while it <numc4:
            if abs(cx0-motormc4[it])<0:
                cx0-=d
                it=-1
            it+=1
        motormc4.append(cx0)


    ran4 = random.uniform(0, 1)
    ran5 = random.uniform(0, 1)
    lxian=l2r-l3l
    lxianp=l1r-l3l
    lxianpx=l2r-l4l
    lspa=l2r-l3l
    if ran4 < ka1* (round(lspa/d)-1-sum(kina))* h:
        cx0=l2r-(round(lspa*random.uniform(0,1)/d)-1)*d
        motor.append([cx0,cx0])
        kina.append(1)
        kinb.append(0)
        ypsl1.append(0)
        dxc.append(0)
        ypsl2.append(0)
        zj0.append(0)
        js0.append(0)
        zj1.append(0)
        js1.append(0)
        zj2.append(0)
        js2.append(0)
        zj3.append(0)
        js3.append(0)
        zj4.append(0)
        js4.append(0)
        zj5.append(0)
        js5.append(0)
        zj6.append(0)
        js6.append(0)
        zj7.append(0)
        js7.append(0)
        pan.append(0)
    if ran5 < ka1* (round(lspa/d)-1-sum(kinb))* h:
        dxc.append(0)
        kina.append(0)
        kinb.append(1)
        cx0=l3l+(round(lspa*random.uniform(0,1)/d)-1)*d
        motor.append([cx0,cx0])
        ypsl1.append(0)
        ypsl2.append(0)
        zj0.append(0)
        js0.append(0)
        zj1.append(0)
        js1.append(0)
        zj2.append(0)
        js2.append(0)
        zj3.append(0)
        js3.append(0)
        zj4.append(0)
        js4.append(0)
        zj5.append(0)
        js5.append(0)
        zj6.append(0)
        js6.append(0)
        zj7.append(0)
        js7.append(0)
        pan.append(0)
    ran4 = random.uniform(0, 1)
    ran5 = random.uniform(0, 1)
    if ran4 < ka2* (round(lxianp/d)-1-sum(kinap))* h:
        cx0=l1r-(round(lxianp*random.uniform(0,1)/d)-1)*d
        motorp.append([cx0,cx0])
        kinap.append(1)
        kinbp.append(0)
        dxcp.append(0)
        ypsl1p.append(0)
        ypsl2p.append(0)
        zj0p.append(0)
        js0p.append(0)
        zj1p.append(0)
        js1p.append(0)
        zj2p.append(0)
        js2p.append(0)
        zj3p.append(0)
        js3p.append(0)
        zj4p.append(0)
        js4p.append(0)
        zj5p.append(0)
        js5p.append(0)
        zj6p.append(0)
        js6p.append(0)
        zj7p.append(0)
        js7p.append(0)
        panp.append(0)   
    if ran5 < ka2* (round(lxianp/d)-1-sum(kinbp))* h:
        cx0=l3l+(round(lxianp*random.uniform(0,1)/d)-1)*d
        motorp.append([cx0,cx0])
        kinap.append(0)
        kinbp.append(1)
        dxcp.append(0)
        ypsl1p.append(0)
        ypsl2p.append(0)               
        panp.append(0)
        zj0p.append(0)
        js0p.append(0)
        zj1p.append(0)
        js1p.append(0)
        zj2p.append(0)
        js2p.append(0)
        zj3p.append(0)
        js3p.append(0)
        zj4p.append(0)
        js4p.append(0)
        zj5p.append(0)
        js5p.append(0)
        zj6p.append(0)
        js6p.append(0)
        zj7p.append(0)
        js7p.append(0)
    ran4 = random.uniform(0, 1)
    ran5 = random.uniform(0, 1)
    if ran4 < ka2* (round(lxianpx/d)-1-sum(kinapx))* h:
        cx0=l4l+(round(lxianpx*random.uniform(0,1)/d)-1)*d
        motorpx.append([cx0,cx0])
        kinapx.append(1)
        kinbpx.append(0)
        dxcpx.append(0)
        ypsl1px.append(0)
        ypsl2px.append(0)
        zj0px.append(0)
        js0px.append(0)
        zj1px.append(0)
        js1px.append(0)
        zj2px.append(0)
        js2px.append(0)
        zj3px.append(0)
        js3px.append(0)
        zj4px.append(0)
        js4px.append(0)
        zj5px.append(0)
        js5px.append(0)
        zj6px.append(0)
        js6px.append(0)
        zj7px.append(0)
        js7px.append(0)
        panpx.append(0)   
    if ran5 < ka2* (round(lxianpx/d)-1-sum(kinbpx))* h:
        cx0=l2r-(round(lxianpx*random.uniform(0,1)/d)-1)*d
        motorpx.append([cx0,cx0])
        kinapx.append(0)
        kinbpx.append(1)
        dxcpx.append(0)
        ypsl1px.append(0)
        ypsl2px.append(0)
        panpx.append(0)
        zj0px.append(0)
        js0px.append(0)
        zj1px.append(0)
        js1px.append(0)
        zj2px.append(0)
        js2px.append(0)
        zj3px.append(0)
        js3px.append(0)
        zj4px.append(0)
        js4px.append(0)
        zj5px.append(0)
        js5px.append(0)
        zj6px.append(0)
        js6px.append(0)
        zj7px.append(0)
        js7px.append(0)
################################################################
    ran4 = random.uniform(0, 1)
    ran5 = random.uniform(0, 1)
    ran4 = random.uniform(0, 1)
    ran5 = random.uniform(0, 1)

    nu=len(kina)
    nu1=len(kinap)
    nu1x=len(kinapx)
    numc1=len(motormc1)
    numc2=len(motormc2)
    numc3=len(motormc3)
    numc4=len(motormc4)
    for i in range(numc3):
        ran4=random.uniform(0, 1)
        if adeps3!=0 and motormc3.index(max(motormc3))==i:
            continue
        if ran4<kdif*h:
            zq=0
            if motormc3[i]+d>=motormc3[motormc3.index(max(motormc3))]-0.1 and  motormc3[i]<motormc3[motormc3.index(max(motormc3))] and adeps3!=0:
                zq=3
            if motormc3[i]+d>l3r:
                zq=3
            for xn in range(numc3):
                if motormc3[i]+d>=motormc3[xn]-0.1 and  motormc3[i]<motormc3[xn]:
                    zq=3
            if zq==0:
                motormc3[i]+=d
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc3[i]-d<l2r:
                zq=3
            for xn in range(numc3):
                if motormc3[i]-d<=motormc3[xn]+0.1 and  motormc3[i]>motormc3[xn]:
                    zq=3
            if zq==0:
                motormc3[i]-=d
    for i in range(numc1):
        if adeps1!=0 and motormc1.index(min(motormc1))==i:
            continue
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc1[i]+d>l3l or motormc1[i]+d>Xkin:
                zq=3
            for xn in range(numc1):
                if motormc1[i]+d>=motormc1[xn]-0.1 and  motormc1[i]<motormc1[xn]:
                    zq=3
            if zq==0:
                motormc1[i]+=d
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc1[i]-d<=motormc1[motormc1.index(min(motormc1))]+0.1 and  motormc1[i]>motormc1[motormc1.index(min(motormc1))] and adeps1!=0:
                zq=3
            if motormc1[i]-d<=l1l:
                zq=3
            for xn in range(numc1):
                if motormc1[i]-d<=motormc1[xn]+0.1 and  motormc1[i]>motormc1[xn]:
                    zq=3
            if zq==0:
                motormc1[i]-=d
    for i in range(numc4):
        ran4=random.uniform(0, 1)
        if adeps4!=0 and motormc4.index(max(motormc4))==i:
            continue
        if ran4<kdif*h:
            zq=0
            if motormc4[i]+d>=motormc4[motormc4.index(max(motormc4))]-0.1 and  motormc4[i]<motormc4[motormc4.index(max(motormc4))] and adeps4!=0:
                zq=3
            if motormc4[i]+d>l4r:
                zq=3
            for xn in range(numc4):
                if motormc4[i]+d>=motormc4[xn]-0.1 and  motormc4[i]<motormc4[xn]:
                    zq=3
            if zq==0:
                motormc4[i]+=d
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc4[i]-d<l2r or motormc4[i]-d<Xkin1:
                zq=3
            for xn in range(numc4):
                if motormc4[i]-d<=motormc4[xn]+0.1 and  motormc4[i]>motormc4[xn]:
                    zq=3
            if zq==0:
                motormc4[i]-=d
    for i in range(numc2):
        if adeps2!=0 and motormc2.index(min(motormc2))==i:
            continue
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc2[i]+d>l3l:
                zq=3
            for xn in range(numc2):
                if motormc2[i]+d>=motormc2[xn]-0.1 and  motormc2[i]<motormc2[xn]:
                    zq=3
            if zq==0:
                motormc2[i]+=d
        ran4=random.uniform(0, 1)
        if ran4<kdif*h:
            zq=0
            if motormc2[i]-d<=motormc2[motormc2.index(min(motormc2))]+0.1 and  motormc2[i]>motormc2[motormc2.index(min(motormc2))] and adeps2!=0:
                zq=3
            if motormc2[i]-d<l2l:
                zq=3
            for xn in range(numc2):
                if motormc2[i]-d<=motormc2[xn]+0.1 and  motormc2[i]>motormc2[xn]:
                    zq=3
            if zq==0:
                motormc2[i]-=d


    if numc3>0: 
        if adeps3!=0 and max(motormc3)+3/2*d<l3r:
            adeps3=0
        if max(motormc3)+3/2*d>=l3r:
            adeps3+=1
            ran4=random.uniform(0, 1)
            if l3r>=Xpo:
                vd=vi
                #stad=stai
            if l3r<Xpo:
                vd=vk
                #stad=stak
            if ran4<vd*h:
                motormc3[motormc3.index(max(motormc3))]-=d
                #print(max(motormc3),l3r,adeps3,'3')
                l3r-=d
                desi3+=d
                pan1=1
                panx=1
            ran4=random.uniform(0, 1)
            if ran4<stad*h:
                motormc3.pop(motormc3.index(max(motormc3)))
                adeps3=0
                numc3-=1
    if numc1>0:
        if adeps1!=0 and min(motormc1)-3/2*d>l1l:
            adeps1=0
        if min(motormc1)-3/2*d<=l1l:
            adeps1+=1
            ran4=random.uniform(0, 1)
            if l1l<=Xpole:
                vd=vi
                #stad=stai
            if l1l>Xpole:
                vd=vk
                #stad=stak
            if ran4<vd*h:
                motormc1[motormc1.index(min(motormc1))]+=d
                pan1=1
                l1l+=d
                desk1+=d
                panx=1
            ran4=random.uniform(0, 1)
            if ran4<stad*h:
                motormc1.pop(motormc1.index(min(motormc1)))
                adeps1=0
                numc1-=1
    if numc4>0: 
        if adeps4!=0 and max(motormc4)+3/2*d<l4r:
            adeps4=0
        if max(motormc4)+3/2*d>=l4r:
            adeps4+=1
            ran4=random.uniform(0, 1)
            if l4r>=Xpo:
                vd=vi
                #stad=stai
            if l4r<Xpo:
                vd=vk
                #stad=stak
            if ran4<vd*h:
                motormc4[motormc4.index(max(motormc4))]-=d
                l4r-=d
                desk4+=d
                pan1=1
                panx=1
            ran4=random.uniform(0, 1)
            if ran4<stad*h:
                motormc4.pop(motormc4.index(max(motormc4)))
                adeps4=0
                numc4-=1
    if numc2>0:
        if adeps2!=0 and min(motormc2)-3/2*d>l2l:
            adeps2=0
        if min(motormc2)-3/2*d<=l2l:
            adeps2+=1
            ran4=random.uniform(0, 1)
            if l2l<=Xpole:
                vd=vi
                #stad=stai
            if l2l>Xpole:
                vd=vk
                #stad=stak
            if ran4<vd*h:
                motormc2[motormc2.index(min(motormc2))]+=d
                #print(min(motormc2),l2l,adeps2,'2')
                pan1=1
                l2l+=d
                desi2+=d
                panx=1
            ran4=random.uniform(0, 1)
            if ran4<stad*h:
                motormc2.pop(motormc2.index(min(motormc2)))
                adeps2=0
                numc2-=1

    ran5 = random.uniform(0, 1)
    vpolkmt=vp0*(1+(Kp3*(Xkin-l1r)/Fp0))/d
    if ran5<abs(vpolkmt)*h:
        if vpolkmt>0:
            l1r+=d
            pos1+=d
        if vpolkmt<0:
            l1r-=d
            pos1-=d
        pan1=1
        panx=1
    ran5 = random.uniform(0, 1)
    #vpolkmt=vp0*(1+(Kp3*(Xkin-l1r)/Fp0))/d
    if ran5<abs(vpolkmt)*h:
        if vpolkmt>0:
            l4l-=d
            pos4-=d
        if vpolkmt<0:
            l4l+=d
            pos4+=d
        pan1=1
        panx=1
    if pan1==1 and panx==1:
        #print('mtd')
        panx=0
        F3=Kp2*(Xpo-l3r)
        F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
        F2=Kp2*(l2l-Xpole)
        F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
        for ton1 in range(nu):
            if kina[ton1]==1 and kinb[ton1]==1:
                dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                if dxc[ton1]>xdmotor:
                    F3+=(dxc[ton1]-xdmotor)*Km
                    F2+=(dxc[ton1]-xdmotor)*Km
                if dxc[ton1]<-xdmotor:
                    F3+=(dxc[ton1]+xdmotor)*Km
                    F2+=(dxc[ton1]+xdmotor)*Km
            if kina[ton1]==0 or kinb[ton1]==0:
                dxc[ton1]=0
        for ton1 in range(nu1):
            if kinap[ton1]==1 and kinbp[ton1]==1:
                dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                if dxcp[ton1]>xdmotor:
                    F3+=(dxcp[ton1]-xdmotor)*Km
                    F1+=(dxcp[ton1]-xdmotor)*Km
                if dxcp[ton1]<-xdmotor:
                    F3+=(dxcp[ton1]+xdmotor)*Km
                    F1+=(dxcp[ton1]+xdmotor)*Km
            if kinap[ton1]==0 or kinbp[ton1]==0:
                dxcp[ton1]=0
        for ton1 in range(nu1x):
            if kinapx[ton1]==1 and kinbpx[ton1]==1:
                dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                if dxcpx[ton1]>xdmotor:
                    F4-=(dxcpx[ton1]-xdmotor)*Km
                    F2-=(dxcpx[ton1]-xdmotor)*Km
                if dxcpx[ton1]<-xdmotor:
                    F4-=(dxcpx[ton1]+xdmotor)*Km
                    F2-=(dxcpx[ton1]+xdmotor)*Km
            if kinapx[ton1]==0 or kinbpx[ton1]==0:
                dxcpx[ton1]=0
        de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
        for xn in range(nu):
            if kina[xn]==1:
                motor[xn][0]+=de1
            if kinb[xn]==1:
                motor[xn][1]+=de2
        for xn in range(nu1):
            if kinap[xn]==1:
                motorp[xn][0]+=de3
            if kinbp[xn]==1:
                motorp[xn][1]+=de2
        for xn in range(nu1x):
            if kinapx[xn]==1:
                motorpx[xn][0]+=de4
            if kinbpx[xn]==1:
                motorpx[xn][1]+=de1
        for xn in range(numc1):
            motormc1[xn]+=de3
        for xn in range(numc2):
            motormc2[xn]+=de1
        for xn in range(numc3):
            motormc3[xn]+=de2
        for xn in range(numc4):
            motormc4[xn]+=de4
        l2l+=de1
        l2r+=de1
        l3l+=de2
        l3r+=de2
        l1l+=de3
        l1r+=de3
        l4l+=de4
        l4r+=de4
        Xpole+=de0
        Xpo+=de7
        Xkin+=de5
        Xkin1+=de6
        Ximt2+=de1
        Ximt3+=de2
        Xkmt1+=de3
        Xkmt4+=de4

    ran4=random.uniform(0, 1)
    if ran4<vpol*h:
        l2r+=d
        if l2r>l3r:
            l2r-=d
    ran4=random.uniform(0, 1)
    if ran4<vpol*h:
        l3l-=d
        if l3l<l2l:
            l3l+=d
    i=0
    while i < numc1:
        if adeps1!=0 and motormc1.index(min(motormc1))==i:
            i=i+1
            continue
        ran4=random.uniform(0, 1)
        if ran4<koff*h:
            del motormc1[i]
            numc1-=1
            continue
        if motormc1[i]<l1l-0.1 or motormc1[i]>l1r+0.1:
            del motormc1[i]
            numc1-=1
            i-=1
        i+=1
    i=0
    while i < numc3:
        if adeps3!=0 and motormc3.index(max(motormc3))==i:
            i=i+1
            continue
        ran4=random.uniform(0, 1)
        if ran4<koff*h:
            del motormc3[i]
            numc3-=1
            continue
        if motormc3[i]<l3l-0.1 or motormc3[i]>l3r+0.1:
            del motormc3[i]
            numc3-=1
            i-=1
        i+=1

    i=0
    while i < numc2:
        if adeps2!=0 and motormc2.index(min(motormc2))==i:
            i=i+1
            continue
        ran4=random.uniform(0, 1)
        if ran4<koff*h:
            del motormc2[i]
            numc2-=1
            continue
        if motormc2[i]<l2l-0.1 or motormc2[i]>l2r+0.1:
            del motormc2[i]
            numc2-=1
            i-=1
        i+=1
    i=0
    while i < numc4:
        if adeps4!=0 and motormc4.index(max(motormc4))==i:
            i=i+1
            continue
        ran4=random.uniform(0, 1)
        if ran4<koff*h:
            del motormc4[i]
            numc4-=1
            continue
        if motormc4[i]<l4l-0.1 or motormc4[i]>l4r+0.1:
            del motormc4[i]
            numc4-=1
            i-=1
        i+=1
#########################################################################              
    nu=len(kina)
    nu1=len(kinap)
    nu1x=len(kinapx)
    monum=sum(kina)
    for xn in range(nu):
        if kina[xn]==1 and kinb[xn]==0:
            monum-=1
    monump=sum(kinap)
    for xn in range(nu1):
        if kinap[xn]==1 and kinbp[xn]==0:
            monump-=1
    monumpx=sum(kinapx)
    for xn in range(nu1x):
        if kinapx[xn]==1 and kinbpx[xn]==0:
            monumpx-=1
    if monum!=0:
        for i in range(0,nu):
            if kina[i]+kinb[i]<2:
                if pan1==1:
                    pan[i]=0
            if kina[i]==1 and kinb[i]==0:
                na=stepnum(0,0,0)
                if na==1:
                    motor[i][0]+=d
                if na==-1:
                    motor[i][0]-=d
            if kina[i]==0 and kinb[i]==1:
                nb=stepnum(0,0,0)
                if nb==1:
                    motor[i][1]-=d
                if nb==-1:
                    motor[i][1]+=d
            if kina[i]==1 and kinb[i]==1:
                if pan1!=0:
                    dx1=(motor[i][0]-motor[i][1])
                    F3=Kp2*(Xpo-l3r)
                    F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                    F2=Kp2*(l2l-Xpole)
                    F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                    pan1=1
                    for ton1 in range(nu):
                        if kina[ton1]==1 and kinb[ton1]==1:
                            dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                            if ton1==i:
                                dxc[i]=dx1+d
                            if dxc[ton1]>xdmotor:
                                F3+=(dxc[ton1]-xdmotor)*Km
                                F2+=(dxc[ton1]-xdmotor)*Km
                            if dxc[ton1]<-xdmotor:
                                F3+=(dxc[ton1]+xdmotor)*Km
                                F2+=(dxc[ton1]+xdmotor)*Km
                        if kina[ton1]==0 or kinb[ton1]==0:
                            dxc[ton1]=0
                    for ton1 in range(nu1):
                        if kinap[ton1]==1 and kinbp[ton1]==1:
                            dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                            if dxcp[ton1]>xdmotor:
                                F3+=(dxcp[ton1]-xdmotor)*Km
                                F1+=(dxcp[ton1]-xdmotor)*Km
                            if dxcp[ton1]<-xdmotor:
                                F3+=(dxcp[ton1]+xdmotor)*Km
                                F1+=(dxcp[ton1]+xdmotor)*Km
                        if kinap[ton1]==0 or kinbp[ton1]==0:
                            dxcp[ton1]=0
                    for ton1 in range(nu1x):
                        if kinapx[ton1]==1 and kinbpx[ton1]==1:
                            dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                            if dxcpx[ton1]>xdmotor:
                                F4-=(dxcpx[ton1]-xdmotor)*Km
                                F2-=(dxcpx[ton1]-xdmotor)*Km
                            if dxcpx[ton1]<-xdmotor:
                                F4-=(dxcpx[ton1]+xdmotor)*Km
                                F2-=(dxcpx[ton1]+xdmotor)*Km
                        if kinapx[ton1]==0 or kinbpx[ton1]==0:
                            dxcpx[ton1]=0
                                            
                    zj0[i],zj1[i],zj2[i],zj3[i],zj4[i],zj5[i],zj6[i],zj7[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxc[i]=dx1-d
                    if dx1+d>xdmotor:
                        F2-=(dx1+d-xdmotor)*Km
                        F3-=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F2-=(dx1+d+xdmotor)*Km
                        F3-=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F2+=(dx1-d-xdmotor)*Km
                        F3+=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F2+=(dx1-d+xdmotor)*Km  
                        F3+=(dx1-d+xdmotor)*Km
                    js0[i],js1[i],js2[i],js3[i],js4[i],js5[i],js6[i],js7[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl0,ypsl1[i],ypsl2[i]=energy(i)

                na=stepnum(ypsl0,ypsl1[i],ypsl2[i])
                nb=stepnum(ypsl0,ypsl1[i],ypsl2[i])
                if motor[i][0]+d>l2r and na==1:
                    na=0
                if motor[i][1]-d<l3l and nb==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na>0.9:
                    zq=0
                    for xn in range(nu):
                        if motor[i][0]+d>=motor[xn][0] and  motor[i][0]<motor[xn][0] and kina[xn]==1:
                            zq=3
                    if zq==0:
                        motor[i][0]+=d
                        pan1=1
                        pan[i]=1
                if na<-0.9:
                    zq=0
                    for xn in range(nu):
                        if motor[i][0]-d<=motor[xn][0] and motor[i][0]>motor[xn][0] and kina[xn]==1:
                            zq=3
                    if zq==0:
                        motor[i][0]-=d
                        pan1=1
                        pan[i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu):
                        if motor[i][1]-d<=motor[xn][1] and  motor[i][1]>motor[xn][1] and kinb[xn]==1:
                            zq=3
                    if zq==0:
                        motor[i][1]-=d
                        pan1=1
                        pan[i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu):
                        if motor[i][1]+d>=motor[xn][1] and  motor[i][1]<motor[xn][1] and kinb[xn]==1:
                            zq=3
                    if zq==0:
                        motor[i][1]+=d
                        pan1=1
                        pan[i]=1
                if pan[i]==1 and zq==0:
                    #print('mz',ypsl0,ypsl1[i],ypsl2[i])
                    if na>0.1 or nb>0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=zj1[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=zj2[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=zj3[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=zj2[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=zj4[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=zj1[i]
                        for xn in range(numc1):
                            motormc1[xn]+=zj3[i]
                        for xn in range(numc2):
                            motormc2[xn]+=zj1[i]
                        for xn in range(numc3):
                            motormc3[xn]+=zj2[i]
                        for xn in range(numc4):
                            motormc4[xn]+=zj4[i]
                        l2l+=zj1[i]
                        l2r+=zj1[i]
                        l3l+=zj2[i]
                        l3r+=zj2[i]
                        l1l+=zj3[i]
                        l1r+=zj3[i]
                        l4l+=zj4[i]
                        l4r+=zj4[i]
                        Xpole+=zj0[i]
                        Xpo+=zj7[i]
                        Xkin+=zj5[i]
                        Xkin1+=zj6[i]
                        Ximt2+=zj1[i]
                        Ximt3+=zj2[i]
                        Xkmt1+=zj3[i]
                        Xkmt4+=zj4[i]
                    if na<-0.1 or nb<-0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=js1[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=js2[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=js3[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=js2[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=js4[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=js1[i]
                        for xn in range(numc1):
                            motormc1[xn]+=js3[i]
                        for xn in range(numc2):
                            motormc2[xn]+=js1[i]
                        for xn in range(numc3):
                            motormc3[xn]+=js2[i]
                        for xn in range(numc4):
                            motormc4[xn]+=js4[i]
                        l2l+=js1[i]
                        l2r+=js1[i]
                        l3l+=js2[i]
                        l3r+=js2[i]
                        l1l+=js3[i]
                        l1r+=js3[i]
                        l4l+=js4[i]
                        l4r+=js4[i]
                        Xpole+=js0[i]
                        Xpo+=js7[i]
                        Xkin+=js5[i]
                        Xkin1+=js6[i]
                        Ximt2+=js1[i]
                        Ximt3+=js2[i]
                        Xkmt1+=js3[i]
                        Xkmt4+=js4[i]
                    ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                    if ckF==1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=de1
                            if kinb[xn]==1:
                                motor[xn][1]+=de2
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=de3
                            if kinbp[xn]==1:
                                motorp[xn][1]+=de2
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=de4
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=de1
                        for xn in range(numc1):
                            motormc1[xn]+=de3
                        for xn in range(numc2):
                            motormc2[xn]+=de1
                        for xn in range(numc3):
                            motormc3[xn]+=de2
                        for xn in range(numc4):
                            motormc4[xn]+=de4
                        l2l+=de1
                        l2r+=de1
                        l3l+=de2
                        l3r+=de2
                        l1l+=de3
                        l1r+=de3
                        l4l+=de4
                        l4r+=de4
                        Xpole+=de0
                        Xpo+=de7
                        Xkin+=de5
                        Xkin1+=de6
                        Ximt2+=de1
                        Ximt3+=de2
                        Xkmt1+=de3
                        Xkmt4+=de4
                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        pan[i]=0


    if monump!=0:
        for i in range(nu1):
            if kinap[i]+kinbp[i]<2:
                panp[i]=0
            if kinap[i]==1 and kinbp[i]==0:
                na=stepnum(0,0,0)
                if na==1:
                    motorp[i][0]+=d
                if na==-1:
                    motorp[i][0]-=d
            if kinap[i]==0 and kinbp[i]==1:
                nb=stepnum(0,0,0)
                if nb==1:
                    motorp[i][1]-=d
                if nb==-1:
                    motorp[i][1]+=d
            if kinap[i]==1 and kinbp[i]==1:
                if pan1!=0:
                    dx1=(motorp[i][0]-motorp[i][1])
                    F3=Kp2*(Xpo-l3r)
                    F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                    F2=Kp2*(l2l-Xpole)
                    F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                    for ton1 in range(nu):
                        if kina[ton1]==1 and kinb[ton1]==1:
                            dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                            if dxc[ton1]>xdmotor:
                                F3+=(dxc[ton1]-xdmotor)*Km
                                F2+=(dxc[ton1]-xdmotor)*Km
                            if dxc[ton1]<-xdmotor:
                                F3+=(dxc[ton1]+xdmotor)*Km
                                F2+=(dxc[ton1]+xdmotor)*Km
                        if kina[ton1]==0 or kinb[ton1]==0:
                            dxc[ton1]=0
                    for ton1 in range(nu1):
                        if kinap[ton1]==1 and kinbp[ton1]==1:
                            dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                            if ton1 ==i:
                                dxcp[i]=dx1+d
                            if dxcp[ton1]>xdmotor:
                                F3+=(dxcp[ton1]-xdmotor)*Km
                                F1+=(dxcp[ton1]-xdmotor)*Km
                            if dxcp[ton1]<-xdmotor:
                                F3+=(dxcp[ton1]+xdmotor)*Km
                                F1+=(dxcp[ton1]+xdmotor)*Km
                        if kinap[ton1]==0 or kinbp[ton1]==0:
                            dxcp[ton1]=0
                    for ton1 in range(nu1x):
                        if kinapx[ton1]==1 and kinbpx[ton1]==1:
                            dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                            if dxcpx[ton1]>xdmotor:
                                F4-=(dxcpx[ton1]-xdmotor)*Km
                                F2-=(dxcpx[ton1]-xdmotor)*Km
                            if dxcpx[ton1]<-xdmotor:
                                F4-=(dxcpx[ton1]+xdmotor)*Km
                                F2-=(dxcpx[ton1]+xdmotor)*Km
                        if kinapx[ton1]==0 or kinbpx[ton1]==0:
                            dxcpx[ton1]=0
                    zj0p[i],zj1p[i],zj2p[i],zj3p[i],zj4p[i],zj5p[i],zj6p[i],zj7p[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxcp[i]=dx1-d
                    if dx1+d>xdmotor:
                        F1-=(dx1+d-xdmotor)*Km
                        F3-=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F1-=(dx1+d+xdmotor)*Km
                        F3-=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F1+=(dx1-d-xdmotor)*Km
                        F3+=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F1+=(dx1-d+xdmotor)*Km   
                        F3+=(dx1-d+xdmotor)*Km
                    js0p[i],js1p[i],js2p[i],js3p[i],js4p[i],js5p[i],js6p[i],js7p[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl0p,ypsl1p[i],ypsl2p[i]=energyp(i)

                na=stepnum(ypsl0p,ypsl1p[i],ypsl2p[i])
                nb=stepnum(ypsl0p,ypsl1p[i],ypsl2p[i])
                if motorp[i][0]+d>l1r and na==1:
                    na=0
                if motorp[i][1]-d<l3l and nb==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na>0.9:
                    zq=0
                    for xn in range(nu1):
                        if motorp[i][0]+d>=motorp[xn][0] and  motorp[i][0]<motorp[xn][0] and kinap[xn]==1:
                            zq=3
                    if zq==0:
                        motorp[i][0]+=d
                        pan1=1
                        panp[i]=1
                if na<-0.9:
                    zq=0
                    for xn in range(nu1):
                        if motorp[i][0]-d<=motorp[xn][0] and motorp[i][0]>motorp[xn][0] and kinap[xn]==1:
                            zq=3
                    if zq==0:
                        motorp[i][0]-=d
                        pan1=1
                        panp[i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu1):
                        if motorp[i][1]-d<=motorp[xn][1] and  motorp[i][1]>motorp[xn][1] and kinbp[xn]==1:
                            zq=3
                    if zq==0:
                        motorp[i][1]-=d
                        pan1=1
                        panp[i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu1):
                        if motorp[i][1]+d>=motorp[xn][1] and  motorp[i][1]<motorp[xn][1] and kinbp[xn]==1:
                            zq=3
                    if zq==0:
                        motorp[i][1]+=d
                        pan1=1
                        panp[i]=1
                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        panp[i]=0
                if panp[i]==1 and zq==0:
                    #print('mpz',ypsl0p,ypsl1p[i],ypsl2p[i],na,nb,motorp[i][0]-motorp[i][1])
                    if na>0.1 or nb>0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=zj1p[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=zj2p[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=zj3p[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=zj2p[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=zj4p[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=zj1p[i]
                        for xn in range(numc1):
                            motormc1[xn]+=zj3p[i]
                        for xn in range(numc2):
                            motormc2[xn]+=zj1p[i]
                        for xn in range(numc3):
                            motormc3[xn]+=zj2p[i]
                        for xn in range(numc4):
                            motormc4[xn]+=zj4p[i]
                        l2l+=zj1p[i]
                        l2r+=zj1p[i]
                        l3l+=zj2p[i]
                        l3r+=zj2p[i]
                        l1l+=zj3p[i]
                        l1r+=zj3p[i]
                        l4l+=zj4p[i]
                        l4r+=zj4p[i]
                        Xpole+=zj0p[i]
                        Xpo+=zj7p[i]
                        Xkin+=zj5p[i]
                        Xkin1+=zj6p[i]
                        Ximt2+=zj1p[i]
                        Ximt3+=zj2p[i]
                        Xkmt1+=zj3p[i]
                        Xkmt4+=zj4p[i]
                    if na<-0.1 or nb<-0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=js1p[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=js2p[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=js3p[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=js2p[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=js4p[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=js1p[i]
                        for xn in range(numc1):
                            motormc1[xn]+=js3p[i]
                        for xn in range(numc2):
                            motormc2[xn]+=js1p[i]
                        for xn in range(numc3):
                            motormc3[xn]+=js2p[i]
                        for xn in range(numc4):
                            motormc4[xn]+=js4p[i]
                        l2l+=js1p[i]
                        l2r+=js1p[i]
                        l3l+=js2p[i]
                        l3r+=js2p[i]
                        l1l+=js3p[i]
                        l1r+=js3p[i]
                        l4l+=js4p[i]
                        l4r+=js4p[i]
                        Xpole+=js0p[i]
                        Xpo+=js7p[i]
                        Xkin+=js5p[i]
                        Xkin1+=js6p[i]
                        Ximt2+=js1p[i]
                        Ximt3+=js2p[i]
                        Xkmt1+=js3p[i]
                        Xkmt4+=js4p[i]


    if monumpx!=0:
        for i in range(nu1x):
            if kinapx[i]+kinbpx[i]<2:
                panpx[i]=0
            if kinapx[i]==1 and kinbpx[i]==0:
                na=stepnum(0,0,0)
                if na==1:
                    motorpx[i][0]-=d
                if na==-1:
                    motorpx[i][0]+=d
            if kinapx[i]==0 and kinbpx[i]==1:
                nb=stepnum(0,0,0)
                if nb==1:
                    motorpx[i][1]+=d
                if nb==-1:
                    motorpx[i][1]-=d
            if kinapx[i]==1 and kinbpx[i]==1:
                if pan1!=0:
                    dx1=(motorpx[i][0]-motorpx[i][1])
                    F3=Kp2*(Xpo-l3r)
                    F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                    F2=Kp2*(l2l-Xpole)
                    F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                    for ton1 in range(nu):
                        if kina[ton1]==1 and kinb[ton1]==1:
                            dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                            if dxc[ton1]>xdmotor:
                                F3+=(dxc[ton1]-xdmotor)*Km
                                F2+=(dxc[ton1]-xdmotor)*Km
                            if dxc[ton1]<-xdmotor:
                                F3+=(dxc[ton1]+xdmotor)*Km
                                F2+=(dxc[ton1]+xdmotor)*Km
                        if kina[ton1]==0 or kinb[ton1]==0:
                            dxc[ton1]=0
                    for ton1 in range(nu1):
                        if kinap[ton1]==1 and kinbp[ton1]==1:
                            dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                            if dxcp[ton1]>xdmotor:
                                F3+=(dxcp[ton1]-xdmotor)*Km
                                F1+=(dxcp[ton1]-xdmotor)*Km
                            if dxcp[ton1]<-xdmotor:
                                F3+=(dxcp[ton1]+xdmotor)*Km
                                F1+=(dxcp[ton1]+xdmotor)*Km
                        if kinap[ton1]==0 or kinbp[ton1]==0:
                            dxcp[ton1]=0
                    for ton1 in range(nu1x):
                        if kinapx[ton1]==1 and kinbpx[ton1]==1:
                            dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                            if ton1==i:
                                dxcpx[ton1]+=d
                            if dxcpx[ton1]>xdmotor:
                                F4-=(dxcpx[ton1]-xdmotor)*Km
                                F2-=(dxcpx[ton1]-xdmotor)*Km
                            if dxcpx[ton1]<-xdmotor:
                                F4-=(dxcpx[ton1]+xdmotor)*Km
                                F2-=(dxcpx[ton1]+xdmotor)*Km
                        if kinapx[ton1]==0 or kinbpx[ton1]==0:
                            dxcpx[ton1]=0
                    zj0px[i],zj1px[i],zj2px[i],zj3px[i],zj4px[i],zj5px[i],zj6px[i],zj7px[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxcpx[i]=dx1-d
                    if dx1+d>xdmotor:
                        F2+=(dx1+d-xdmotor)*Km
                        F4+=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F2+=(dx1+d+xdmotor)*Km
                        F4+=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F2-=(dx1-d-xdmotor)*Km
                        F4-=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F2-=(dx1-d+xdmotor)*Km   
                        F4-=(dx1-d+xdmotor)*Km
                    js0px[i],js1px[i],js2px[i],js3px[i],js4px[i],js5px[i],js6px[i],js7px[i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl0px,ypsl1px[i],ypsl2px[i]=energypx(i)

                na=stepnum(ypsl0px,ypsl2px[i],ypsl1px[i])
                nb=stepnum(ypsl0px,ypsl2px[i],ypsl1px[i]) 
                if motorpx[i][1]+d>l2r and nb==1:
                    na=0
                if motorpx[i][0]-d<l4l and na==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na<-0.9:
                    zq=0
                    for xn in range(nu1x):
                        if motorpx[i][0]+d>=motorpx[xn][0] and  motorpx[i][0]<motorpx[xn][0] and kinapx[xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[i][0]+=d
                        pan1=1
                        panpx[i]=1

                if na>0.9:
                    zq=0
                    for xn in range(nu1x):
                        if motorpx[i][0]-d<=motorpx[xn][0] and motorpx[i][0]>motorpx[xn][0] and kinapx[xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[i][0]-=d
                        pan1=1
                        panpx[i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu1x):
                        if motorpx[i][1]-d<=motorpx[xn][1] and  motorpx[i][1]>motorpx[xn][1] and kinbpx[xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[i][1]-=d
                        pan1=1
                        panpx[i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu1x):
                        if motorpx[i][1]+d>=motorpx[xn][1] and  motorpx[i][1]<motorpx[xn][1] and kinbpx[xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[i][1]+=d
                        pan1=1
                        panpx[i]=1
                if panpx[i]==1 and zq==0:
                    #print('mpxz',ypsl0px,ypsl2px[i],ypsl1px[i],dxcpx[i])
                    if na<-0.1 or nb<-0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=zj1px[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=zj2px[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=zj3px[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=zj2px[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=zj4px[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=zj1px[i]
                        for xn in range(numc1):
                            motormc1[xn]+=zj3px[i]
                        for xn in range(numc2):
                            motormc2[xn]+=zj1px[i]
                        for xn in range(numc3):
                            motormc3[xn]+=zj2px[i]
                        for xn in range(numc4):
                            motormc4[xn]+=zj4px[i]
                        l2l+=zj1px[i]
                        l2r+=zj1px[i]
                        l3l+=zj2px[i]
                        l3r+=zj2px[i]
                        l1l+=zj3px[i]
                        l1r+=zj3px[i]
                        l4l+=zj4px[i]
                        l4r+=zj4px[i]
                        Xpole+=zj0px[i]
                        Xpo+=zj7px[i]
                        Xkin+=zj5px[i]
                        Xkin1+=zj6px[i]
                        Ximt2+=zj1px[i]
                        Ximt3+=zj2px[i]
                        Xkmt1+=zj3px[i]
                        Xkmt4+=zj4px[i]
                    if na>0.1 or nb>0.1:
                        for xn in range(nu):
                            if kina[xn]==1:
                                motor[xn][0]+=js1px[i]
                            if kinb[xn]==1:
                                motor[xn][1]+=js2px[i]
                        for xn in range(nu1):
                            if kinap[xn]==1:
                                motorp[xn][0]+=js3px[i]
                            if kinbp[xn]==1:
                                motorp[xn][1]+=js2px[i]
                        for xn in range(nu1x):
                            if kinapx[xn]==1:
                                motorpx[xn][0]+=js4px[i]
                            if kinbpx[xn]==1:
                                motorpx[xn][1]+=js1px[i]
                        for xn in range(numc1):
                            motormc1[xn]+=js3px[i]
                        for xn in range(numc2):
                            motormc2[xn]+=js1px[i]
                        for xn in range(numc3):
                            motormc3[xn]+=js2px[i]
                        for xn in range(numc4):
                            motormc4[xn]+=js4px[i]
                        l2l+=js1px[i]
                        l2r+=js1px[i]
                        l3l+=js2px[i]
                        l3r+=js2px[i]
                        l1l+=js3px[i]
                        l1r+=js3px[i]
                        l4l+=js4px[i]
                        l4r+=js4px[i]
                        Xpole+=js0px[i]
                        Xpo+=js7px[i]
                        Xkin+=js5px[i]
                        Xkin1+=js6px[i]
                        Ximt2+=js1px[i]
                        Ximt3+=js2px[i]
                        Xkmt1+=js3px[i]
                        Xkmt4+=js4px[i]

                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        panpx[i]=0          
    if sum(pan)+sum(panp)+sum(panpx)==0:
        pan1=0
###############################################       
##########################################################  
    tuogeng=0
    for xbb in range(nu):
        xpa=0
        if kina[xbb]==0 and kinb[xbb]==1:
            xpa=1
            kinb[xbb]=dissociation(0,0,0,0)
            if motor[xbb][1]>=l2l and motor[xbb][1]<=l2r:
                kina[xbb]=rebind()
            if kinb[xbb]==1 and kina[xbb]==1:
                lian1=0
                for xn in range(nu):
                    if kina[xn]==1 and xn!=xbb and abs(l2r-round((l2r-motor[xbb][1])/d)*d-motor[xn][0])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l2r+round((motor[xbb][1]-l2r)/d)*d>=motor[xbb][1]:
                        la=(round((motor[xbb][1]-l2r)/d)-1)*d
                    if l2r+round((motor[xbb][1]-l2r)/d)*d<motor[xbb][1]:
                        la=(round((motor[xbb][1]-l2r)/d)+1)*d
                    for xn in range(nu):
                        if kina[xn]==1 and xn!=xbb and abs(l2r+la-motor[xn][0])<1:
                            lian1=2
                            break
                if lian1==0:
                    motor[xbb][0]=l2r+round((motor[xbb][1]-l2r)/d)*d
                if lian1==1:
                    motor[xbb][0]=l2r+la
                if lian1==2:
                    kina[xbb]=0
                pan1=1
                monum+=1
        if xpa==0 and kina[xbb]==1 and kinb[xbb]==0:
            xpa=1
            kina[xbb]=dissociation(0,0,0,0)
            if motor[xbb][0]>=l3l and motor[xbb][0]<=l3r:
                kinb[xbb]=rebind()
            if kinb[xbb]==1 and kina[xbb]==1:
                lian1=0
                for xn in range(nu):
                    if kinb[xn]==1 and xn!=xbb and abs(l3l+round((motor[xbb][0]-l3l)/d)*d-motor[xn][1])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l3l+round((motor[xbb][0]-l3l)/d)*d>=motor[xbb][0]:
                        la=(round((motor[xbb][0]-l3l)/d)-1)*d
                    if l3l+round((motor[xbb][0]-l3l)/d)*d<motor[xbb][0]:
                        la=(round((motor[xbb][0]-l3l)/d)+1)*d
                    for xn in range(nu):
                        if kinb[xn]==1 and xn!=xbb and abs(l3l+la-motor[xn][1])<1:
                            lian1=2
                            break
                if lian1==0:
                    motor[xbb][1]=l3l+round((motor[xbb][0]-l3l)/d)*d
                if lian1==1:
                    motor[xbb][1]=l3l+la
                if lian1==2:
                    kinb[xbb]=0
                pan1=1
                monum+=1
        if xpa==0 and kina[xbb]==1 and kinb[xbb]==1:
            F=0
            if motor[xbb][0]-motor[xbb][1]>xdmotor:
                F=Km*(motor[xbb][0]-motor[xbb][1]-xdmotor)
            if motor[xbb][0]-motor[xbb][1]<-xdmotor:
                F=Km*(motor[xbb][0]-motor[xbb][1]+xdmotor)
            kina[xbb]=dissociation(F,ypsl0px,ypsl1[xbb],ypsl2[xbb])
            kinb[xbb]=dissociation(F,ypsl0px,ypsl1[xbb],ypsl2[xbb])
            if motor[xbb][0]>l2r or motor[xbb][0]<l2l:
                kina[xbb]=0
            if motor[xbb][1]<l3l or motor[xbb][1]>l3r:
                kinb[xbb]=0
            if kina[xbb]==0 or kinb[xbb]==0:
                monum-=1
                pan1=1
                F3=Kp2*(Xpo-l3r)
                F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                F2=Kp2*(l2l-Xpole)
                F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                for ton1 in range(nu):
                    if kina[ton1]==1 and kinb[ton1]==1:
                        dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                        if dxc[ton1]>xdmotor:
                            F3+=(dxc[ton1]-xdmotor)*Km
                            F2+=(dxc[ton1]-xdmotor)*Km
                        if dxc[ton1]<-xdmotor:
                            F3+=(dxc[ton1]+xdmotor)*Km
                            F2+=(dxc[ton1]+xdmotor)*Km
                    if kina[ton1]==0 or kinb[ton1]==0:
                        dxc[ton1]=0
                for ton1 in range(nu1):
                    if kinap[ton1]==1 and kinbp[ton1]==1:
                        dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                        if dxcp[ton1]>xdmotor:
                            F3+=(dxcp[ton1]-xdmotor)*Km
                            F1+=(dxcp[ton1]-xdmotor)*Km
                        if dxcp[ton1]<-xdmotor:
                            F3+=(dxcp[ton1]+xdmotor)*Km
                            F1+=(dxcp[ton1]+xdmotor)*Km
                    if kinap[ton1]==0 or kinbp[ton1]==0:
                        dxcp[ton1]=0
                for ton1 in range(nu1x):
                    if kinapx[ton1]==1 and kinbpx[ton1]==1:
                        dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                        if dxcpx[ton1]>xdmotor:
                            F4-=(dxcpx[ton1]-xdmotor)*Km
                            F2-=(dxcpx[ton1]-xdmotor)*Km
                        if dxcpx[ton1]<-xdmotor:
                            F4-=(dxcpx[ton1]+xdmotor)*Km
                            F2-=(dxcpx[ton1]+xdmotor)*Km
                    if kinapx[ton1]==0 or kinbpx[ton1]==0:
                        dxcpx[ton1]=0
                de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                
                for xn in range(nu):
                    if kina[xn]==1:
                        motor[xn][0]+=de1
                    if kinb[xn]==1:
                        motor[xn][1]+=de2
                for xn in range(nu1):
                    if kinap[xn]==1:
                        motorp[xn][0]+=de3
                    if kinbp[xn]==1:
                        motorp[xn][1]+=de2
                for xn in range(nu1x):
                    if kinapx[xn]==1:
                        motorpx[xn][0]+=de4
                    if kinbpx[xn]==1:
                        motorpx[xn][1]+=de1
                for xn in range(numc1):
                    motormc1[xn]+=de3
                for xn in range(numc2):
                    motormc2[xn]+=de1
                for xn in range(numc3):
                    motormc3[xn]+=de2
                for xn in range(numc4):
                    motormc4[xn]+=de4
                l2l+=de1
                l2r+=de1
                l3l+=de2
                l3r+=de2
                l1l+=de3
                l1r+=de3
                l4l+=de4
                l4r+=de4
                Xpole+=de0
                Xpo+=de7
                Xkin+=de5
                Xkin1+=de6
                Ximt2+=de1
                Ximt3+=de2
                Xkmt1+=de3
                Xkmt4+=de4
                ypsl0px,ypsl1[xbb],ypsl2[xbb]=energy(xbb)
                #print('mtuo')
        if kina[xbb]==0 and kinb[xbb]==0:
            tuogeng=1
######################################################
    if tuogeng==1:
        for xn in range(nu):
            if kina[xn]==0 and kinb[xn]==0:
                tuogeng=0
                del motor[xn]
                del kina[xn]
                del kinb[xn]
                del dxc[xn]
                del ypsl1[xn]
                del ypsl2[xn]
                del zj7[xn]
                del js7[xn]
                del zj6[xn]
                del js6[xn]
                del zj5[xn]
                del js5[xn]
                del zj4[xn]
                del js4[xn]
                del zj3[xn]
                del js3[xn]
                del zj2[xn]
                del js2[xn]
                del zj1[xn]
                del js1[xn]
                del zj0[xn]
                del js0[xn]
                del pan[xn]
                nu-=1
                break 
##################################################################
    tuogeng=0
    nu=len(kina)
    for xbb in range(nu1):
        xpa=0 
        if kinap[xbb]==0 and kinbp[xbb]==1:
            xpa=1
            kinbp[xbb]=dissociation(0,0,0,0)
            if motorp[xbb][1]>=l3l and motorp[xbb][1]<=l1r:
                kinap[xbb]=rebind()
            if kinbp[xbb]==1 and kinap[xbb]==1:
                lian1=0
                for xn in range(nu1):
                    if kinap[xn]==1 and xn!=xbb and abs(l1l-round((l1l-motorp[xbb][1])/d)*d-motorp[xn][0])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l1l+round((motorp[xbb][1]-l1l)/d)*d>=motorp[xbb][1]:
                        la=(round((motorp[xbb][1]-l1l)/d)-1)*d
                    if l1l+round((motorp[xbb][1]-l1l)/d)*d<motorp[xbb][1]:
                        la=(round((motorp[xbb][1]-l1l)/d)+1)*d
                    for xn in range(nu1):
                        if kinap[xn]==1 and xn!=xbb and abs(l1l+la-motorp[xn][0])<1:
                            lian1=2
                            break
                if lian1==0:
                    motorp[xbb][0]=l1l+round((motorp[xbb][1]-l1l)/d)*d
                if lian1==1:
                    motorp[xbb][0]=l1l+la
                if lian1==2:
                    kinap[xbb]=0
                pan1=1
                monump+=1
        if xpa==0 and kinap[xbb]==1 and kinbp[xbb]==0:
            xpa=1
            kinap[xbb]=dissociation(0,0,0,0)
            if motorp[xbb][0]>=l3l and motorp[xbb][0]<=l1r:
                kinbp[xbb]=rebind()
            if kinbp[xbb]==1 and kinap[xbb]==1:
                lian1=0
                for xn in range(nu1):
                    if kinbp[xn]==1 and xn!=xbb and abs(l3l+round((motorp[xbb][0]-l3l)/d)*d-motorp[xn][1])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l3l+round((motorp[xbb][0]-l3l)/d)*d>=motorp[xbb][0]:
                        la=(round((motorp[xbb][0]-l3l)/d)-1)*d
                    if l3l+round((motorp[xbb][0]-l3l)/d)*d<motorp[xbb][0]:
                        la=(round((motorp[xbb][0]-l3l)/d)+1)*d
                    for xn in range(nu1):
                        if kinbp[xn]==1 and xn!=xbb and abs(l3l+la-motorp[xn][1])<1:
                            lian1=2
                            break
                if lian1==0:
                    motorp[xbb][1]=l3l+round((motorp[xbb][0]-l3l)/d)*d
                if lian1==1:
                    motorp[xbb][1]=l3l+la
                if lian1==2:
                    kinbp[xbb]=0
                pan1=1
                monump+=1
        if xpa==0 and  kinap[xbb]==1 and kinbp[xbb]==1:
            F=0
            if motorp[xbb][0]-motorp[xbb][1]>xdmotor:
                F=Km*(motorp[xbb][0]-motorp[xbb][1]-xdmotor)
            if motorp[xbb][0]-motorp[xbb][1]<-xdmotor:
                F=Km*(motorp[xbb][0]-motorp[xbb][1]+xdmotor)
            kinap[xbb]=dissociation(F,ypsl0px,ypsl1p[xbb],ypsl2p[xbb])
            kinbp[xbb]=dissociation(F,ypsl0px,ypsl1p[xbb],ypsl2p[xbb])
            if motorp[xbb][0]>l1r or motorp[xbb][0]<l1l:
                kinap[xbb]=0
            if motorp[xbb][1]<l3l or motorp[xbb][1]>l3r:
                kinbp[xbb]=0
            if kinap[xbb]==0 or kinbp[xbb]==0:
                monump-=1
                pan1=1
                F3=Kp2*(Xpo-l3r)
                F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                F2=Kp2*(l2l-Xpole)
                F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                for ton1 in range(nu):
                    if kina[ton1]==1 and kinb[ton1]==1:
                        dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                        if dxc[ton1]>xdmotor:
                            F3+=(dxc[ton1]-xdmotor)*Km
                            F2+=(dxc[ton1]-xdmotor)*Km
                        if dxc[ton1]<-xdmotor:
                            F3+=(dxc[ton1]+xdmotor)*Km
                            F2+=(dxc[ton1]+xdmotor)*Km
                    if kina[ton1]==0 or kinb[ton1]==0:
                        dxc[ton1]=0
                for ton1 in range(nu1):
                    if kinap[ton1]==1 and kinbp[ton1]==1:
                        dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                        if dxcp[ton1]>xdmotor:
                            F3+=(dxcp[ton1]-xdmotor)*Km
                            F1+=(dxcp[ton1]-xdmotor)*Km
                        if dxcp[ton1]<-xdmotor:
                            F3+=(dxcp[ton1]+xdmotor)*Km
                            F1+=(dxcp[ton1]+xdmotor)*Km
                    if kinap[ton1]==0 or kinbp[ton1]==0:
                        dxcp[ton1]=0
                for ton1 in range(nu1x):
                    if kinapx[ton1]==1 and kinbpx[ton1]==1:
                        dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                        if dxcpx[ton1]>xdmotor:
                            F4-=(dxcpx[ton1]-xdmotor)*Km
                            F2-=(dxcpx[ton1]-xdmotor)*Km
                        if dxcpx[ton1]<-xdmotor:
                            F4-=(dxcpx[ton1]+xdmotor)*Km
                            F2-=(dxcpx[ton1]+xdmotor)*Km
                    if kinapx[ton1]==0 or kinbpx[ton1]==0:
                        dxcpx[ton1]=0
                de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                for xn in range(nu):
                    if kina[xn]==1:
                        motor[xn][0]+=de1
                    if kinb[xn]==1:
                        motor[xn][1]+=de2
                for xn in range(nu1):
                    if kinap[xn]==1:
                        motorp[xn][0]+=de3
                    if kinbp[xn]==1:
                        motorp[xn][1]+=de2
                for xn in range(nu1x):
                    if kinapx[xn]==1:
                        motorpx[xn][0]+=de4
                    if kinbpx[xn]==1:
                        motorpx[xn][1]+=de1
                for xn in range(numc1):
                    motormc1[xn]+=de3
                for xn in range(numc2):
                    motormc2[xn]+=de1
                for xn in range(numc3):
                    motormc3[xn]+=de2
                for xn in range(numc4):
                    motormc4[xn]+=de4
                l2l+=de1
                l2r+=de1
                l3l+=de2
                l3r+=de2
                l1l+=de3
                l1r+=de3
                l4l+=de4
                l4r+=de4
                Xpole+=de0
                Xpo+=de7
                Xkin+=de5
                Xkin1+=de6
                Ximt2+=de1
                Ximt3+=de2
                Xkmt1+=de3
                Xkmt4+=de4
                ypsl0px,ypsl1p[xbb],ypsl2p[xbb]=energyp(xbb)  
                #print('mptuo')
        if kinap[xbb]==0 and kinbp[xbb]==0:
            tuogeng=1
    if tuogeng==1:
        for xn in range(nu1):
            if kinap[xn]==0 and kinbp[xn]==0:
                tuogeng=0
                del motorp[xn]
                del kinap[xn]
                del kinbp[xn]
                del dxcp[xn]
                del ypsl1p[xn]
                del ypsl2p[xn]
                del zj7p[xn]
                del js7p[xn]
                del zj6p[xn]
                del js6p[xn]
                del zj5p[xn]
                del js5p[xn]
                del zj4p[xn]
                del js4p[xn]
                del zj3p[xn]
                del js3p[xn]
                del zj2p[xn]
                del js2p[xn]
                del zj1p[xn]
                del js1p[xn]
                del zj0p[xn]
                del js0p[xn]
                del panp[xn]
                nu1-=1
                break
    tuogeng=0
    for xbb in range(nu1x):
        xpa=0
        if kinapx[xbb]==0 and kinbpx[xbb]==1:
            xpa=1
            kinbpx[xbb]=dissociation(0,0,0,0)
            if motorpx[xbb][1]>=l4l and motorpx[xbb][1]<=l2r:
                kinapx[xbb]=rebind()

            if kinbpx[xbb]==1 and kinapx[xbb]==1:
                lian1=0
                #print('111')
                for xn in range(nu1x):
                    if kinapx[xn]==1 and xn!=xbb and abs(l4l-round((l4l-motorpx[xbb][1])/d)*d-motorpx[xn][0])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l4l+round((motorpx[xbb][1]-l4l)/d)*d>=motorpx[xbb][1]:
                        la=(round((motorpx[xbb][1]-l4l)/d)-1)*d
                    if l4l+round((motorpx[xbb][1]-l4l)/d)*d<motorpx[xbb][1]:
                        la=(round((motorpx[xbb][1]-l4l)/d)+1)*d
                    for xn in range(nu1x):
                        if kinapx[xn]==1 and xn!=xbb and abs(l4l+la-motorpx[xn][0])<1:
                            lian1=2
                            break
                if lian1==0:
                    motorpx[xbb][0]=l4l+round((motorpx[xbb][1]-l4l)/d)*d
                if lian1==1:
                    motorpx[xbb][0]=l4l+la
                if lian1==2:
                    kinapx[xbb]=0
                    #print('222')
                pan1=1
                monumpx+=1
        if xpa==0 and kinapx[xbb]==1 and kinbpx[xbb]==0:
            xpa=1
            kinapx[xbb]=dissociation(0,0,0,0)
            if motorpx[xbb][0]>=l4l and motorpx[xbb][0]<=l2r:
                kinbpx[xbb]=rebind()
            if kinbpx[xbb]==1 and kinapx[xbb]==1:
                lian1=0
                for xn in range(nu1x):
                    if kinbpx[xn]==1 and xn!=xbb and abs(l2l+round((motorpx[xbb][0]-l2l)/d)*d-motorpx[xn][1])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l2l+round((motorpx[xbb][0]-l2l)/d)*d>=motorpx[xbb][0]:
                        la=(round((motorpx[xbb][0]-l2l)/d)-1)*d
                    if l2l+round((motorpx[xbb][0]-l2l)/d)*d<motorpx[xbb][0]:
                        la=(round((motorpx[xbb][0]-l2l)/d)+1)*d
                    for xn in range(nu1x):
                        if kinbpx[xn]==1 and xn!=xbb and abs(l2l+la-motorpx[xn][1])<1:
                            lian1=2
                            break
                if lian1==0:
                    motorpx[xbb][1]=l2l+round((motorpx[xbb][0]-l2l)/d)*d
                if lian1==1:
                    motorpx[xbb][1]=l2l+la
                if lian1==2:
                    kinbpx[xbb]=0
                pan1=1
                monumpx+=1
        if xpa==0 and kinapx[xbb]==1 and kinbpx[xbb]==1:
            F=0
            if motorpx[xbb][0]-motorpx[xbb][1]>xdmotor:
                F=-Km*(motorpx[xbb][0]-motorpx[xbb][1]-xdmotor)
            if motorpx[xbb][0]-motorpx[xbb][1]<-xdmotor:
                F=-Km*(motorpx[xbb][0]-motorpx[xbb][1]+xdmotor)
            kinapx[xbb]=dissociation(F,ypsl0px,ypsl2px[xbb],ypsl1px[xbb])
            kinbpx[xbb]=dissociation(F,ypsl0px,ypsl2px[xbb],ypsl1px[xbb])
            if motorpx[xbb][0]>l4r or motorpx[xbb][0]<l4l:
                kinapx[xbb]=0
            if motorpx[xbb][1]<l2l or motorpx[xbb][1]>l2r:
                kinbpx[xbb]=0
            if kinapx[xbb]==0 or kinbpx[xbb]==0:
                monumpx-=1
                pan1=1
                F3=Kp2*(Xpo-l3r)
                F1=Kp1*(l1l-Xpole)+Kp3*(l1r-Xkin)
                F2=Kp2*(l2l-Xpole)
                F4=Kp1*(Xpo-l4r)+Kp3*(Xkin1-l4l)
                for ton1 in range(nu):
                    if kina[ton1]==1 and kinb[ton1]==1:
                        dxc[ton1]=motor[ton1][0]-motor[ton1][1]
                        if dxc[ton1]>xdmotor:
                            F3+=(dxc[ton1]-xdmotor)*Km
                            F2+=(dxc[ton1]-xdmotor)*Km
                        if dxc[ton1]<-xdmotor:
                            F3+=(dxc[ton1]+xdmotor)*Km
                            F2+=(dxc[ton1]+xdmotor)*Km
                    if kina[ton1]==0 or kinb[ton1]==0:
                        dxc[ton1]=0
                for ton1 in range(nu1):
                    if kinap[ton1]==1 and kinbp[ton1]==1:
                        dxcp[ton1]=motorp[ton1][0]-motorp[ton1][1]
                        if dxcp[ton1]>xdmotor:
                            F3+=(dxcp[ton1]-xdmotor)*Km
                            F1+=(dxcp[ton1]-xdmotor)*Km
                        if dxcp[ton1]<-xdmotor:
                            F3+=(dxcp[ton1]+xdmotor)*Km
                            F1+=(dxcp[ton1]+xdmotor)*Km
                    if kinap[ton1]==0 or kinbp[ton1]==0:
                        dxcp[ton1]=0
                for ton1 in range(nu1x):
                    if kinapx[ton1]==1 and kinbpx[ton1]==1:
                        dxcpx[ton1]=motorpx[ton1][0]-motorpx[ton1][1]
                        if dxcpx[ton1]>xdmotor:
                            F4-=(dxcpx[ton1]-xdmotor)*Km
                            F2-=(dxcpx[ton1]-xdmotor)*Km
                        if dxcpx[ton1]<-xdmotor:
                            F4-=(dxcpx[ton1]+xdmotor)*Km
                            F2-=(dxcpx[ton1]+xdmotor)*Km
                    if kinapx[ton1]==0 or kinbpx[ton1]==0:
                        dxcpx[ton1]=0
                de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                for xn in range(nu):
                    if kina[xn]==1:
                        motor[xn][0]+=de1
                    if kinb[xn]==1:
                        motor[xn][1]+=de2
                for xn in range(nu1):
                    if kinap[xn]==1:
                        motorp[xn][0]+=de3
                    if kinbp[xn]==1:
                        motorp[xn][1]+=de2
                for xn in range(nu1x):
                    if kinapx[xn]==1:
                        motorpx[xn][0]+=de4
                    if kinbpx[xn]==1:
                        motorpx[xn][1]+=de1
                for xn in range(numc1):
                    motormc1[xn]+=de3
                for xn in range(numc2):
                    motormc2[xn]+=de1
                for xn in range(numc3):
                    motormc3[xn]+=de2
                for xn in range(numc4):
                    motormc4[xn]+=de4
                l2l+=de1
                l2r+=de1
                l3l+=de2
                l3r+=de2
                l1l+=de3
                l1r+=de3
                l4l+=de4
                l4r+=de4
                Xpole+=de0
                Xpo+=de7
                Xkin+=de5
                Xkin1+=de6
                Ximt2+=de1
                Ximt3+=de2
                Xkmt1+=de3
                Xkmt4+=de4
                ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                if ckF==1:
                    for xn in range(nu):
                        if kina[xn]==1:
                            motor[xn][0]+=de1
                        if kinb[xn]==1:
                            motor[xn][1]+=de2
                    for xn in range(nu1):
                        if kinap[xn]==1:
                            motorp[xn][0]+=de3
                        if kinbp[xn]==1:
                            motorp[xn][1]+=de2
                    for xn in range(nu1x):
                        if kinapx[xn]==1:
                            motorpx[xn][0]+=de4
                        if kinbpx[xn]==1:
                            motorpx[xn][1]+=de1
                    for xn in range(numc1):
                        motormc1[xn]+=de3
                    for xn in range(numc2):
                        motormc2[xn]+=de1
                    for xn in range(numc3):
                        motormc3[xn]+=de2
                    for xn in range(numc4):
                        motormc4[xn]+=de4
                    l2l+=de1
                    l2r+=de1
                    l3l+=de2
                    l3r+=de2
                    l1l+=de3
                    l1r+=de3
                    l4l+=de4
                    l4r+=de4
                    Xpole+=de0
                    Xpo+=de7
                    Xkin+=de5
                    Xkin1+=de6
                    Ximt2+=de1
                    Ximt3+=de2
                    Xkmt1+=de3
                    Xkmt4+=de4
                ypsl0px,ypsl1px[xbb],ypsl2px[xbb]=energypx(xbb)
                
        if kinapx[xbb]==0 and kinbpx[xbb]==0:
            tuogeng=1
    if tuogeng==1:
        for xn in range(nu1x):
            if kinapx[xn]==0 and kinbpx[xn]==0:
                tuogeng=0
                del motorpx[xn]
                del kinapx[xn]
                del kinbpx[xn]
                del dxcpx[xn]
                del ypsl1px[xn]
                del ypsl2px[xn]
                del zj7px[xn]
                del js7px[xn]
                del zj6px[xn]
                del js6px[xn]
                del zj5px[xn]
                del js5px[xn]
                del zj4px[xn]
                del js4px[xn]
                del zj3px[xn]
                del js3px[xn]
                del zj2px[xn]
                del js2px[xn]
                del zj1px[xn]
                del js1px[xn]
                del zj0px[xn]
                del js0px[xn]
                del panpx[xn]
                nu1x=len(kinapx)
                break

    if lxian<0 or Xpole>0:
        break
    if a%1000==0:
        time1.append(a*h)
        MTX1.append(lxian)
        MTX2.append(round(Xpo-Xpole,2))
        MTzn6.append(monum)
        MTzn7.append(monump)
        MTzn8.append(monumpx)
        MTzn9.append(numc1)
        MTzn10.append(numc2)
        MTzn11.append(numc3)
        MTzn12.append(numc4)
        MTzn13.append(round(Xkin,2))
        MTzn14.append(round(Xkin1,2))
        MTzn15.append(round((l2r+l3l)/2,2))
        MTk1.append(round((xkmta1-Xkmt1),2))
        MTi1.append(round((ximta2-Ximt2),2))
        MTzn1.append(round((desi2-desia2),2))
        MTzn2.append(round((desk1-deska1),2))
        MTzn3.append(round((pos1-posa1),2)) 
               
        MTzn4.append(-Kp2*(l2l-Xpole))
        MTzn5.append(Kp3*(Xkin-l1r))

        desia2=desi2
        desia3=desi3
        deska1=desk1
        deska4=desk4
        posa1=pos1
        posa4=pos4
        ximta2=Ximt2
        ximta3=Ximt3
        xkmta1=Xkmt1
        xkmta4=Xkmt4
    if a%100000==0:
        dataframe = pd.DataFrame({'atime':time1,'lover':MTX1,'lpole':MTX2,'pmtk':MTk1,'pmti':MTi1,'rdi2':MTzn1,'rdk1':MTzn2,'rp1':MTzn3,'r4F':MTzn4,'r5F':MTzn5,'n1':MTzn6,'n2':MTzn7,'n3':MTzn8,'n4c':MTzn9,'n5c':MTzn10,'n6c':MTzn11,'n7c':MTzn12,'xkinc':MTzn13,'xkin1c':MTzn14,'xpcp':MTzn14})
        dataframe.to_csv('./2'+str(kon)+str(koff)+str(nant)+str(npar)+str(vi)+str(vk)+str(stai)+str(pbei)+str(Fp0)+str(Kp1)+str(vp1)+'.csv',index=False,sep=',')






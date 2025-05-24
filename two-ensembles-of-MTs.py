#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
    if ypsl0-ypsl1>1000:
        ypsl1=ypsl0-1000
    if ypsl0-ypsl1<-1000:
        ypsl1=ypsl0+1000
    if ypsl0-ypsl2>1000:
        ypsl2=ypsl0-1000
    if ypsl0-ypsl2<-1000:
        ypsl2=ypsl0+1000
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
    p=1
    delta=1
    kr=4.6
    kw0=5
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
    de=[0,[0 for i in range(pair)],[0 for i in range(pair)],[0 for i in range(pair)],[0 for i in range(pair)],0,0,0]
    num1=[0 for i in range(pair)]
    num2=[0 for i in range(pair)]
    nump1=[0 for i in range(pair)]
    nump2=[0 for i in range(pair)]
    numpx1=[0 for i in range(pair)]
    numpx2=[0 for i in range(pair)]
    fen0=[[] for i in range(pair)]
    fen1=[[] for i in range(pair)]
    fen2=[[] for i in range(pair)]
    fenp0=[[] for i in range(pair)]
    fenp1=[[] for i in range(pair)]
    fenp2=[[] for i in range(pair)]
    fenpx0=[[] for i in range(pair)]
    fenpx1=[[] for i in range(pair)]
    fenpx2=[[] for i in range(pair)]
    jisuan=[[] for i in range(pair)]
    jisuanp=[[] for i in range(pair)]
    jisuanpx=[[] for i in range(pair)]
    ji=[[] for i in range(pair)]
    F0=0
    F7=0
    F5=0
    F6=0
    Fcc=0
    for x2 in range(pair):
        F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
        F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
        F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
        if abs(F1[x2])>0.00001 or abs(F2[x2])>0.00001 or abs(F3[x2])>0.00001 or abs(F4[x2])>0.00001:
            Fcc=1
    if abs(F0)>0.0001 or abs(F5)>0.00001 or abs(F6)>0.00001 or abs(F7)>0.00001:
        Fcc=1
    if Fcc==0:
        return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]
    #print(F1,F2,F3,F4,F0,F5,F6,F7)
    for x2 in range(pair):
        nu[x2]=len(kina[x2])
        nu1[x2]=len(kinap[x2])
        nu1x[x2]=len(kinapx[x2])
        for ton1 in range(nu[x2]):
            if kina[x2][ton1]==1 and kinb[x2][ton1]==1:
                if dxc[x2][ton1]>xdmotor:
                    num1[x2]+=1
                if dxc[x2][ton1]<-xdmotor:
                    num2[x2]+=1
        for ton1 in range(nu1[x2]):
            if kinap[x2][ton1]==1 and kinbp[x2][ton1]==1:
                if dxcp[x2][ton1]>xdmotor:
                    nump1[x2]+=1
                if dxcp[x2][ton1]<-xdmotor:
                    nump2[x2]+=1
        for ton1 in range(nu1x[x2]):
            if kinapx[x2][ton1]==1 and kinbpx[x2][ton1]==1:
                if dxcpx[x2][ton1]>xdmotor:
                    numpx1[x2]+=1
                if dxcpx[x2][ton1]<-xdmotor:
                    numpx2[x2]+=1
        jisuan[x2]=sorted(dxc[x2])
        for ton1 in range(nu[x2]-monum[x2]):
            jisuan[x2].remove(0)
        jisuanp[x2]=sorted(dxcp[x2])
        for ton1 in range(nu1[x2]-monump[x2]):
            jisuanp[x2].remove(0)
        jisuanpx[x2]=sorted(dxcpx[x2])
        for ton1 in range(nu1x[x2]-monumpx[x2]):
            jisuanpx[x2].remove(0)
        cnn0=0
        zq=0

        m1=[monum[x2]-num1[x2]]
        m0=[num2[x2]]
        m2=[0]
        if monum[x2]>0:
            ttt=jisuan[x2][0]
            for tt1 in range(monum[x2]):
                if abs(jisuan[x2][tt1]-ttt)>0.0001:
                    if jisuan[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuan[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuan[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuan[x2][tt1]
            m2.append(num2[x2])
            m0.append(monum[x2]-num1[x2])
            m1.append(monum[x2])
        fen2[x2]=list(set(m2))
        fen1[x2]=list(set(m1))
        fen0[x2]=list(set(m0))
        fen2[x2].sort(key=m2.index)
        fen1[x2].sort(key=m1.index)
        fen0[x2].sort(key=m0.index)

        m1=[monump[x2]-nump1[x2]]
        m0=[nump2[x2]]
        m2=[0]
        if monump[x2]>0:
            ttt=jisuanp[x2][0]
            for tt1 in range(monump[x2]):
                if abs(jisuanp[x2][tt1]-ttt)>0.0001:
                    if jisuanp[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuanp[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuanp[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuanp[x2][tt1]
            m2.append(nump2[x2])
            m0.append(monump[x2]-nump1[x2])
            m1.append(monump[x2])
        fenp2[x2]=list(set(m2))
        fenp1[x2]=list(set(m1))
        fenp0[x2]=list(set(m0))
        fenp2[x2].sort(key=m2.index)
        fenp1[x2].sort(key=m1.index)
        fenp0[x2].sort(key=m0.index)

        m1=[monumpx[x2]-numpx1[x2]]
        m0=[numpx2[x2]]
        m2=[0]
        if monumpx[x2]>0:
            ttt=jisuanpx[x2][0]
            for tt1 in range(monumpx[x2]):
                if abs(jisuanpx[x2][tt1]-ttt)>0.0001:
                    if jisuanpx[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuanpx[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuanpx[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuanpx[x2][tt1]
            m2.append(numpx2[x2])
            m0.append(monumpx[x2]-numpx1[x2])
            m1.append(monumpx[x2])
        fenpx2[x2]=list(set(m2))
        fenpx1[x2]=list(set(m1))
        fenpx0[x2]=list(set(m0))
        fenpx2[x2].sort(key=m2.index)
        fenpx1[x2].sort(key=m1.index)
        fenpx0[x2].sort(key=m0.index)    
    #########################################算dzk,dz0
    
    dzk=[[[0 for xxn2 in range(len(fen0[0]))] for xxn1 in range(2)]]
    dz0=[[[0 for xxn2 in range(len(fen2[0]))],[0 for xxn1 in range(len(fen1[0]))]]]
    dzkp=[[[0 for xxn2 in range(len(fenp0[0]))] for xxn1 in range(2)]]
    dz0p=[[[0 for xxn2 in range(len(fenp1[0]))],[0 for xxn1 in range(len(fenp2[0]))]]]
    dzkpx=[[[0 for xxn2 in range(len(fenpx0[0]))] for xxn1 in range(2)]]
    dz0px=[[[0 for xxn2 in range(len(fenpx1[0]))],[0 for xxn1 in range(len(fenpx2[0]))]]]
    for x3 in range(1,pair):
        dzk.append([[0 for xxn2 in range(len(fen0[x3]))] for xxn1 in range(2)])
        dz0.append([[0 for xxn2 in range(len(fen2[x3]))],[0 for xxn1 in range(len(fen1[x3]))]])
        dzkp.append([[0 for xxn2 in range(len(fenp0[x3]))] for xxn1 in range(2)])
        dz0p.append([[0 for xxn2 in range(len(fenp1[x3]))],[0 for xxn1 in range(len(fenp2[x3]))]])
        dzkpx.append([[0 for xxn2 in range(len(fenpx0[x3]))] for xxn1 in range(2)])
        dz0px.append([[0 for xxn2 in range(len(fenpx1[x3]))],[0 for xxn1 in range(len(fenpx2[x3]))]])
    for x2 in range(pair):
        ta1=-1
        for xck in fen0[x2]:
            ta1+=1
            cnnk=monum[x2]-num1[x2]-fen0[x2][-ta1-1]
            chek=monum[x2]-num1[x2]-1#dxc增，ds1>0
            for xxxn in range(cnnk):
                dzk[x2][0][ta1]+=jisuan[x2][chek-xxxn]-xdmotor
        ta1=-1
        for xck in fen0[x2]:
            ta1+=1
            cnnk=xck-num2[x2]
            chek=num2[x2]
            for xxxn in range(cnnk):
                dzk[x2][1][ta1]+=jisuan[x2][chek+xxxn]+xdmotor 
        ta2=-1
        for xckp in fenp0[x2]:
            ta2+=1
            cnnkp=xckp-nump2[x2]
            chekp=nump2[x2]
            for xxxn in range(cnnkp):
                dzkp[x2][0][ta2]+=jisuanp[x2][chekp+xxxn]+xdmotor
        ta2=-1
        for xckp in fenp0[x2]:
            ta2+=1
            cnnkp=monump[x2]-nump1[x2]-fenp0[x2][-ta2-1]
            chekp=monump[x2]-nump1[x2]-1 
            for xxxn in range(cnnkp):
                dzkp[x2][1][ta2]+=jisuanp[x2][chekp-xxxn]-xdmotor 
        ta2=-1
        for xckpx in fenpx0[x2]:
            ta2+=1
            cnnkpx=xckpx-numpx2[x2]
            chekpx=numpx2[x2]
            for xxxn in range(cnnkpx):
                dzkpx[x2][0][ta2]+=jisuanpx[x2][chekpx+xxxn]+xdmotor
        ta2=-1
        for xckpx in fenpx0[x2]:
            ta2+=1
            cnnkpx=monumpx[x2]-numpx1[x2]-fenpx0[x2][-ta2-1]
            chekpx=monumpx[x2]-numpx1[x2]-1 
            for xxxn in range(cnnkpx):
                dzkpx[x2][1][ta2]+=jisuanpx[x2][chekpx-xxxn]-xdmotor 
        ta3=-1
        for xc0 in fen2[x2]:
            ta3+=1
            cnn0=num2[x2]-fen2[x2][-ta3-1]
            che0=num2[x2]-1
            for xxxn in range(cnn0):
                dz0[x2][0][ta3]+=jisuan[x2][che0-xxxn]+xdmotor
        ta3=-1
        for xc0 in fen1[x2]:    
            ta3+=1
            cnn0=xc0-monum[x2]+num1[x2]
            che0=monum[x2]-num1[x2]
            for xxxn in range(cnn0):
                dz0[x2][1][ta3]+=jisuan[x2][che0+xxxn]-xdmotor
        ta4=-1
        for xc0p in fenp1[x2]:
            ta4+=1
            cnn0p=xc0p-monump[x2]+nump1[x2]
            che0p=monump[x2]-nump1[x2]
            for xxxn in range(cnn0p):#dxcp减，ds3>0
                dz0p[x2][0][ta4]+=jisuanp[x2][che0p+xxxn]-xdmotor
        ta4=-1
        for xc0p in fenp2[x2]:
            ta4+=1
            cnn0p=nump2[x2]-fenp2[x2][-ta4-1]
            che0p=nump2[x2]-1
            for xxxn in range(cnn0p):
                dz0p[x2][1][ta4]+=jisuanp[x2][che0p-xxxn]+xdmotor
        ta4=-1
        for xc0px in fenpx1[x2]:
            ta4+=1
            cnn0px=xc0px-monumpx[x2]+numpx1[x2]
            che0px=monumpx[x2]-numpx1[x2]
            for xxxn in range(cnn0px):
                dz0px[x2][0][ta4]+=jisuanpx[x2][che0px+xxxn]-xdmotor
        ta4=-1
        for xc0px in fenpx2[x2]:
            ta4+=1
            cnn0px=numpx2[x2]-fenpx2[x2][-ta4-1]
            che0px=numpx2[x2]-1
            for xxxn in range(cnn0px):
                dz0px[x2][1][ta4]+=jisuanpx[x2][che0px-xxxn]+xdmotor

    #print(round(F0),round(F1,2),round(F2,2),round(F3,2),round(F4,2),round(F5,2),round(F6,2),round(F7,2))
    for x2 in range(pair):
        for xx1 in range(2):
            ta1=-1
            for xck in fen0[x2]:
                ta1+=1
                if xx1==0: 
                    if ta1<len(dzk[x2][0]):
                        cnnk=monum[x2]-num1[x2]-fen0[x2][-ta1-1]
                    if ta1>=len(dzk[x2][0]):
                        break
                if xx1==1:
                    if ta1<len(dzk[x2][1]):
                        cnnk=xck-num2[x2]
                    if ta1>=len(dzk[x2][1]):
                        break
                for xx2 in range(2):
                    ta2=-1
                    for xckp in fenp0[x2]:
                        ta2+=1
                        if xx2==0:
                            if ta2<len(dzkp[x2][0]):
                                cnnkp=xckp-nump2[x2]
                            if ta2>=len(dzkp[x2][0]):
                                break
                        if xx2==1:
                            if ta2<len(dzkp[x2][1]):
                                cnnkp=monump[x2]-nump1[x2]-fenp0[x2][-ta2-1]
                            if ta2>=len(dzkp[x2][1]):
                                break
                        for xx3 in range(2):
                            ta5=-1
                            for xckpx in fenpx0[x2]:    
                                ta5+=1
                                if xx3==0:
                                    if ta5<len(dzkpx[x2][0]):
                                        cnnkpx=xckpx-numpx2[x2]
                                    if ta5>=len(dzkpx[x2][0]):
                                        break
                                if xx3==1:
                                    if ta5<len(dzkpx[x2][1]):
                                        cnnkpx=monumpx[x2]-numpx1[x2]-fenpx0[x2][-ta5-1]
                                    if ta5>=len(dzkpx[x2][1]):
                                        break
                                if xx1==0:
                                    xa3=fen2[x2]
                                if xx1==1:
                                    xa3=fen1[x2]
                                ta3=-1
                                for xc0 in xa3:
                                    #print(xqq,xc0,len(dz0[xqq][1]),len(dz0[xqq][0]),'3')
                                    ta3+=1  
                                    if xx1==0:
                                        if ta3<len(dz0[x2][0]):
                                            cnn0=num2[x2]-fen2[x2][-ta3-1]
                                        if ta3>=len(dz0[x2][0]):
                                            break
                                    if xx1==1:
                                        if ta3<len(dz0[x2][1]):
                                            cnn0=xc0-monum[x2]+num1[x2]
                                        if ta3>=len(dz0[x2][1]):
                                            break
                                    if xx2==0:
                                        xa4=fenp1[x2]
                                    if xx2==1:
                                        xa4=fenp2[x2]
                                    ta4=-1
                                    for xc0p in xa4:
                                        ta4+=1
                                        if xx2==0:
                                            if ta4<len(dz0p[x2][0]):
                                                cnn0p=xc0p-monump[x2]+nump1[x2]
                                            if ta4>=len(dz0p[x2][0]):
                                                break
                                        if xx2==1:
                                            if ta4<len(dz0p[x2][1]):
                                                cnn0p=nump2[x2]-fenp2[x2][-ta4-1]    
                                        if xx3==0:
                                            xa5=fenpx1[x2]
                                        if xx3==1:
                                            xa5=fenpx2[x2]
                                        ta6=-1
                                        for xc0px in xa5:
                                            ta6+=1
                                            if xx3==0:
                                                if ta6<len(dz0px[x2][0]):
                                                    cnn0px=xc0px-monumpx[x2]+numpx1[x2]
                                                if ta6>=len(dz0px[x2][0]):
                                                    break
                                            if xx3==1:
                                                if ta6<len(dz0px[x2][1]):
                                                    cnn0px=numpx2[x2]-fenpx2[x2][-ta6-1] 
                                                if ta6>len(dz0px[x2][1]):
                                                    break
                                            if cnn0+cnnk==0 and xx1!=0:
                                                continue
                                            if cnn0p+cnnkp==0 and xx2!=0:
                                                continue
                                            if cnn0px+cnnkpx==0 and xx3!=0:
                                                continue    
                                            if cnn0+cnnk!=0:
                                                if xx1==0  and F2[x2]>0.00001 and abs(F3[x2])<0.00001:
                                                    continue
                                                if xx1==1 and F2[x2]<-0.00001 and abs(F3[x2])<0.00001:
                                                    continue
                                                if xx1==0 and F2[x2]>0.00001 and F3[x2]>0.00001:
                                                    continue
                                                if xx1==1  and F2[x2]<-0.00001 and F2[x2]<-0.00001:
                                                    continue
                                                if xx1==0  and F3[x2]>0.00001 and abs(F2[x2])<0.00001:
                                                    continue
                                                if xx1==1 and F3[x2]<-0.00001 and abs(F2[x2])<0.00001:
                                                    continue
                                            if cnn0p+cnnkp!=0:
                                                if xx2==1 and F1[x2]>0.00001 and F3[x2]>=-0.00001:
                                                    continue
                                                if xx2==0 and F1[x2]<-0.00001 and F3[x2]<=0.00001:
                                                    continue
                                                if xx2==1 and F3[x2]>0.00001 and abs(F1[x2])<0.00001:
                                                    continue
                                                if xx2==0 and F3[x2]<-0.00001 and abs(F1[x2])<0.00001:
                                                    continue
                                                if xx2==1 and F1[x2]>0.00001 and abs(F3[x2])<0.00001:
                                                    continue
                                                if xx2==0 and F1[x2]<-0.00001 and abs(F3[x2])<0.00001:
                                                    continue
                                            if cnn0px+cnnkpx!=0:
                                                if xx3==0 and F4[x2]>0.00001 and F2[x2]>=-0.00001:
                                                    continue
                                                if xx3==1 and F4[x2]<-0.00001 and F2[x2]<=0.00001:
                                                    continue
                                                if xx3==0 and F2[x2]>0.00001 and abs(F4[x2])<0.00001:
                                                    continue
                                                if xx3==1 and F2[x2]<-0.00001 and abs(F4[x2])<0.00001:
                                                    continue
                                                if xx3==0 and F4[x2]>0.00001 and abs(F2[x2])<0.00001:
                                                    continue
                                                if xx3==1 and F4[x2]<-0.00001 and abs(F2[x2])<0.00001:
                                                    continue
                                            ji[x2].append([cnn0,cnnk,cnn0p,cnnkp,cnn0px,cnnkpx,xx1,xx2,xx3,ta1,ta2,ta3,ta4,ta5,ta6])
    a00=round(2*(-Kp1-Kp2),3);a10=Kp2;a30=Kp1;a80=Kp2;a100=Kp1
    a21=Kp2;a41=Kp1;a71=round((-Kp1-Kp2)*2,3);a91=Kp2;a111a=Kp1
    a32=Kp3;a52=round((-Kp3-Kp4)*2,3);a62=round(Kp4*2,3);a102=Kp3
    a43=Kp3;a53=round(Kp4*2,3);a63=round((-Kp3-Kp4 )*2,3);a113=Kp3
    zq=[0 for i in range(pair)]
    cxx=[0 for i in range(pair)]
    for cx0 in range(len(ji[0])):
        cxx[0]=cx0
        if abs(sum(zq)-pair)<0.01:
            break
        for cx1 in range(len(ji[1])):
            cxx[1]=cx1
            if abs(sum(zq)-pair)<0.01:
                break
             
            cnn0=ji[0][cxx[0]][0]
            cnnk=ji[0][cxx[0]][1]
            cnn0p=ji[0][cxx[0]][2]
            cnnkp=ji[0][cxx[0]][3]
            cnn0px=ji[0][cxx[0]][4]
            cnnkpx=ji[0][cxx[0]][5]
            xx1=ji[0][cxx[0]][6]
            xx2=ji[0][cxx[0]][7]
            xx3=ji[0][cxx[0]][8]
            ta1=ji[0][cxx[0]][9]
            ta2=ji[0][cxx[0]][10]
            ta3=ji[0][cxx[0]][11]
            ta4=ji[0][cxx[0]][12]
            ta5=ji[0][cxx[0]][13]
            ta6=ji[0][cxx[0]][14]
            a04=-Kp1;a24=round(-(nump1[0]+nump2[0]+cnnkp-cnn0p)*Km,3);a34=round((nump1[0]+nump2[0]+cnnkp-cnn0p)*Km+Kp1+Kp3+Kpn5*nunum1[0],3);
            a54=-Kp3;a104=round(-Kpn5*nunum1[0],3)
            a05=-Kp2;a15=round(((numpx1[0]+numpx2[0]+cnnkpx-cnn0px)+(num1[0]+num2[0]+cnnk-cnn0))*Km+Kp2+Kpn5*nunum[0],3);
            a25=round(-Km*(num1[0]+num2[0]+cnnk-cnn0),3);a45=round(-(numpx1[0]+numpx2[0]+cnnkpx-cnn0px)*Km,3);a85=round(-Kpn5*nunum[0],3)
            a16=round((num1[0]+num2[0]+cnnk-cnn0)*Km,3);a26=round(-Km*(nump1[0]+nump2[0]+cnnkp-cnn0p+num1[0]+num2[0]+cnnk-cnn0)-Kp2-Kpn5*nunum[1],3)
            a36=round(Km*(nump1[0]+nump2[0]+cnnkp-cnn0p),3);a76=Kp2;a96=round(Kpn5*nunum[1],3)
            a17=round(Km*(numpx1[0]+numpx2[0]+cnnkpx-cnn0px),3);a47=round(-Km*(numpx1[0]+numpx2[0]+cnnkpx-cnn0px)-Kp1-Kp3-Kpn5*nunum1[1],3)
            a67=Kp3;a77=Kp1;a117=round(Kpn5*nunum1[1],3)
            b4=-Km*(dzkp[0][xx2][ta2]-dz0p[0][xx2][ta4])-F1[0]
            b5=-Km*(-dzkpx[0][xx3][ta5]+dz0px[0][xx3][ta6]+dzk[0][xx1][ta1]-dz0[0][xx1][ta3])-F2[0]
            b6=-Km*(dzkp[0][xx2][ta2]-dz0p[0][xx2][ta4]+dzk[0][xx1][ta1]-dz0[0][xx1][ta3])-F3[0]
            b7=Km*(dzkpx[0][xx3][ta5]-dz0px[0][xx3][ta6])-F4[0]
            cnn0=ji[1][cxx[1]][0]
            cnnk=ji[1][cxx[1]][1]
            cnn0p=ji[1][cxx[1]][2]
            cnnkp=ji[1][cxx[1]][3]
            cnn0px=ji[1][cxx[1]][4]
            cnnkpx=ji[1][cxx[1]][5]
            xx1=ji[1][cxx[1]][6]
            xx2=ji[1][cxx[1]][7]
            xx3=ji[1][cxx[1]][8]
            ta1=ji[1][cxx[1]][9]
            ta2=ji[1][cxx[1]][10]
            ta3=ji[1][cxx[1]][11]
            ta4=ji[1][cxx[1]][12]
            ta5=ji[1][cxx[1]][13]
            ta6=ji[1][cxx[1]][14]
            a08=-Kp1;a28=round(-(nump1[1]+nump2[1]+cnnkp-cnn0p)*Km,3);a38=round((nump1[1]+nump2[1]+cnnkp-cnn0p)*Km+Kp1+Kp3+Kpn5*nunum1[0],3)
            a58=-Kp3;a038=-round(Kpn5*nunum1[0],3)
            a09=-Kp2;a19=round(((numpx1[1]+numpx2[1]+cnnkpx-cnn0px)+(num1[1]+num2[1]+cnnk-cnn0))*Km+Kp2+Kpn5*nunum[0],3);a019=-round(Kpn5*nunum[0],3)
            a29=round(-Km*(num1[1]+num2[1]+cnnk-cnn0),3);a49=round(-(numpx1[1]+numpx2[1]+cnnkpx-cnn0px)*Km,3);
            a110=round((num1[1]+num2[1]+cnnk-cnn0)*Km,3);a210=round(-Km*(nump1[1]+nump2[1]+cnnkp-cnn0p+num1[1]+num2[1]+cnnk-cnn0)-Kp2-Kpn5*nunum[1],3)
            a310=round(Km*(nump1[1]+nump2[1]+cnnkp-cnn0p),3);a710=Kp2;a0210=round(Kpn5*nunum[1],3)
            a111=round(Km*(numpx1[1]+numpx2[1]+cnnkpx-cnn0px),3);a411=round(-Km*(numpx1[1]+numpx2[1]+cnnkpx-cnn0px)-Kp1-Kp3-Kpn5*nunum1[1],3)
            a611=Kp3;a711=Kp1;a0411=round(Kpn5*nunum1[1],3)
            b8=-Km*(dzkp[1][xx2][ta2]-dz0p[1][xx2][ta4])-F1[1]
            b9=-Km*(-dzkpx[1][xx3][ta5]+dz0px[1][xx3][ta6]+dzk[1][xx1][ta1]-dz0[1][xx1][ta3])-F2[1]
            b10=-Km*(dzkp[1][xx2][ta2]-dz0p[1][xx2][ta4]+dzk[1][xx1][ta1]-dz0[1][xx1][ta3])-F3[1]
            b11=Km*(dzkpx[1][xx3][ta5]-dz0px[1][xx3][ta6])-F4[1] 
            b=np.array([round(-F0,14),round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14),round(b8,14),round(b9,14),round(b10,14),round(b11,14)])
            try:
                A=np.array([[a00,a10,0,a30,0,0,0,0,a80,0,a100,0],[0,0,a21,0,a41,0,0,a71,0,a91,0,a111a],
                                [0,0,0,a32,0,a52,a62,0,0,0,a102,0],[0,0,0,0,a43,a53,a63,0,0,0,0,a113],
                                [a04,0,a24,a34,0,a54,0,0,0,0,a104,0],[a05,a15,a25,0,a45,0,0,0,a85,0,0,0],
                                [0,a16,a26,a36,0,0,0,a76,0,a96,0,0],[0,a17,0,0,a47,0,a67,a77,0,0,0,a117],
                                [a08,0,0,a038,0,a58,0,0,0,a28,a38,0],[a09,a019,0,0,0,0,0,0,a19,a29,0,a49],
                                [0,0,a0210,0,0,0,0,a710,a110,a210,a310,0],[0,0,0,0,a0411,0,a611,a711,a111,0,0,a411]])
                #b=np.array([-F0,-F7,-F5,-F6,b4,b5,b6,b7])

                mtui=np.linalg.solve(A,b)

            except:
                mtui = gaussian_elimination(A,b)
            tuzx=(mtui[0]+mtui[7])/2
            for xn in range(12):
                mtui[xn]-=tuzx

            de=[mtui[0],[mtui[1],mtui[8]],[mtui[2],mtui[9]],[mtui[3],mtui[10]],[mtui[4],mtui[11]],mtui[5],mtui[6],mtui[7]]
            #print(de,'xx')
            #print(b)
            chezq=0
            for x3 in range(pair):
                check0=-1################
                checkk=-1###############
                checkp0=-1###############
                checkpk=-1###############
                checkpx0=-1###############
                checkpxk=-1##############
                cnn0=ji[x3][cxx[x3]][0]
                cnnk=ji[x3][cxx[x3]][1]
                cnn0p=ji[x3][cxx[x3]][2]
                cnnkp=ji[x3][cxx[x3]][3]
                cnn0px=ji[x3][cxx[x3]][4]
                cnnkpx=ji[x3][cxx[x3]][5]
                xx1=ji[x3][cxx[x3]][6]
                xx2=ji[x3][cxx[x3]][7]
                xx3=ji[x3][cxx[x3]][8]
                #print(chezq)
                if cnn0+cnnk==0:
                    for ton1 in range(monum[x3]):
                        if abs(jisuan[x3][ton1])<=xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuan[x3][ton1])>xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0+cnnk!=0:    
                    if xx1==0 and de[1][x3]-de[2][x3]>=0:
                        check0=0###############
                        checkk=0###############
                        for ton1 in range(monum[x3]):
                            if jisuan[x3][ton1]<-xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                                check0+=1
                            if abs(jisuan[x3][ton1])<=xdmotor and jisuan[x3][ton1]+de[1][x3]-de[2][x3]>xdmotor:
                                checkk+=1
                    if xx1==1 and de[1][x3]-de[2][x3]<=0:
                        check0=0###############
                        checkk=0###############
                        for ton1 in range(monum[x3]):
                            if abs(jisuan[x3][ton1])<=xdmotor and jisuan[x3][ton1]+de[1][x3]-de[2][x3]<-xdmotor:
                                checkk+=1
                            if jisuan[x3][ton1]>xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                                check0+=1
                    if check0!=cnn0 or checkk!=cnnk:
                        chezq=-1
                        break
                #print(chezq)
                if chezq==-1:
                    #print(x3)
                    break
                if cnn0p+cnnkp==0:
                    for ton1 in range(monump[x3]):
                        if abs(jisuanp[x3][ton1])<=xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuanp[x3][ton1])>xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0p+cnnkp!=0:
                    if xx2==1 and de[3][x3]-de[2][x3]>=0:
                        checkp0=0###############
                        checkpk=0###############
                        for ton1 in range(monump[x3]):
                            if jisuanp[x3][ton1]<-xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                                checkp0+=1
                            if abs(jisuanp[x3][ton1])<=xdmotor and jisuanp[x3][ton1]+de[3][x3]-de[2][x3]>xdmotor:
                                checkpk+=1
                    if xx2==0 and de[3][x3]-de[2][x3]<=0:#
                        checkp0=0
                        checkpk=0
                        for ton1 in range(monump[x3]):
                            if abs(jisuanp[x3][ton1])<=xdmotor and jisuanp[x3][ton1]+de[3][x3]-de[2][x3]<-xdmotor:
                                checkpk+=1
                            if jisuanp[x3][ton1]>xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                                checkp0+=1
                    if checkp0!=cnn0p or checkpk!=cnnkp:
                        chezq=-1
                        break
                #print(chezq)
                if chezq==-1:
                    #print(x3,'p')
                    break
                if cnn0px+cnnkpx==0:
                    for ton1 in range(monumpx[x3]):
                        if abs(jisuanpx[x3][ton1])<=xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuanpx[x3][ton1])>xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0px+cnnkpx!=0:
                    if xx3==1 and de[4][x3]-de[1][x3]>=0:
                        checkpx0=0################
                        checkpxk=0################
                        for ton1 in range(monumpx[x3]):
                            if jisuanpx[x3][ton1]<-xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                                checkpx0+=1
                            if abs(jisuanpx[x3][ton1])<=xdmotor and jisuanpx[x3][ton1]+de[4][x3]-de[1][x3]>xdmotor:
                                checkpxk+=1
                    if xx3==0 and de[4][x3]-de[1][x3]<=0:
                        checkpx0=0###############
                        checkpxk=0###############
                        for ton1 in range(monumpx[x3]):
                            if abs(jisuanpx[x3][ton1])<=xdmotor and jisuanpx[x3][ton1]+de[4][x3]-de[1][x3]<-xdmotor:
                                checkpxk+=1
                            if jisuanpx[x3][ton1]>xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                                checkpx0+=1
                    if checkpx0!=cnn0px or checkpxk!=cnnkpx:
                        chezq=-1
                        break
                #print(chezq)
                if chezq==-1:
                    #print(x3,'px')
                    break
            if chezq==-1:
                continue
            #print(jisuan,jisuanp,jisuanpx)
            #print(ji[0][cxx[0]],ji[1][cxx[1]])
            #print(ji)
            return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]
    #print(jisuan,jisuanp,jisuanpx)
    dpn=0
    if sum(zq)!=pair:
        for xqq in range(pair):
            if abs(F1[xqq])>0.001 or abs(F2[xqq])>0.001 or abs(F3[xqq])>0.001 or abs(F4[xqq])>0.001 :
                dpn=xqq
                break
    #print(dpn,num1[dpn]+numpx2[dpn]+num2[dpn]+numpx1[dpn]+1,nump2[dpn]+nump1[dpn]+1)
    tq=1
    if abs(F2[dpn])>0.0001 and abs(F3[dpn])<0.0001:
        tq=0
    p1=[0 for xxx in range(pair)]
    p2=[0 for xxx in range(pair)]
    p3=[0 for xxx in range(pair)]
    p4=[0 for xxx in range(pair)]
    for xxx in range(pair):
        p1[xxx]=F1[xxx]
        p2[xxx]=F2[xxx]
        p3[xxx]=F3[xxx]
        p4[xxx]=F4[xxx]
    if tq==0:
        lastui1=-1000
        for xn in range(num1[dpn]+numpx2[dpn]+num2[dpn]+numpx1[dpn]+1):#MT2
            if xn < num1[dpn]:
                tui1=-(jisuan[dpn][-num1[dpn]+xn]-xdmotor+0.1)
            if xn >= num1[dpn] and xn <num1[dpn]+numpx2[dpn]:
                tui1=(jisuanpx[dpn][numpx2[dpn]-1-xn+num1[dpn]]+xdmotor-0.1)
            if xn >= num1[dpn]+numpx2[dpn] and xn <num1[dpn]+numpx2[dpn]+num2[dpn]:
                tui1=-(jisuan[dpn][num2[dpn]-1-xn+num1[dpn]+numpx2[dpn]]+xdmotor-0.1)
            if xn >= num1[dpn]+numpx2[dpn]+num2[dpn] and xn <num1[dpn]+numpx2[dpn]+num2[dpn]+numpx1[dpn]:
                tui1=(jisuanpx[dpn][monumpx[dpn]-numpx1[dpn]+xn-num1[dpn]-numpx2[dpn]-num2[dpn]]-xdmotor+0.1)
            if xn==num1[dpn]+numpx2[dpn]+num2[dpn]+numpx1[dpn]:
                tui1=0
            if abs(tui1-lastui1)<0.0001:
                continue
            lastui1=tui1
            lastui2=-1000
            for xn1 in range(nump2[dpn]+nump1[dpn]+1):
                if xn1 <nump1[dpn]:
                    tui2=(-jisuanp[dpn][-nump1[dpn]+xn1]+xdmotor-0.1)
                if xn1 >= nump1[dpn] and xn1 <nump1[dpn]+nump2[dpn]:
                    tui2=-(jisuanp[dpn][nump2[dpn]-1-xn1+nump1[dpn]]+xdmotor-0.1)
                if xn1==nump1[dpn]+nump2[dpn]:
                    tui2=0
                if abs(tui2-lastui2)<0.0001 or abs(tui1)+abs(tui2)==0:
                    continue
                lastui2=tui2#####MT3-move-tui1,MT4-move-tui2
                F0=0
                F7=0
                F5=0
                F6=0
                for x2 in range(pair):
                    if x2!=dpn:
                        F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                        F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
                        F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(monum[x2]):
                            if jisuan[x2][ton1]>xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                            if jisuan[x2][ton1]<-xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                        for ton1 in range(monump[x2]):
                            if jisuanp[x2][ton1]>xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                            if jisuanp[x2][ton1]<-xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                        for ton1 in range(monumpx[x2]):
                            if jisuanpx[x2][ton1]>xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                            if jisuanpx[x2][ton1]<-xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                    if x2==dpn:
                        F0+=Kp1*(l1l[x2]+tui2-Xpole)+Kp2*(l2l[x2]+tui1-Xpole)
                        F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
                        F5+=Kp3*(l1r[x2]+tui2-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]+tui2-Xpole)+Kp3*(l1r[x2]+tui2-Xkin)
                        F2[x2]=Kp2*(l2l[x2]+tui1-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(monum[x2]):
                            if jisuan[x2][ton1]+tui1>xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]+tui1-xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]+tui1-xdmotor)*Km
                            if jisuan[x2][ton1]+tui1<-xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]+tui1+xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]+tui1+xdmotor)*Km
                        for ton1 in range(monump[x2]):
                            if jisuanp[x2][ton1]+tui2>xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]+tui2-xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]+tui2-xdmotor)*Km
                            if jisuanp[x2][ton1]+tui2<-xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]+tui2+xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]+tui2+xdmotor)*Km
                        for ton1 in range(monumpx[x2]):
                            if jisuanpx[x2][ton1]-tui1>xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]-tui1-xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]-tui1-xdmotor)*Km
                            if jisuanpx[x2][ton1]-tui1<-xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]-tui1+xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]-tui1+xdmotor)*Km
                for ton1 in range(Nma[0]):
                    if dpn==0:
                        kua0=tui1;kua1=0
                    if dpn==1:
                        kua0=0;kua1=tui1 
                    if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                        F2[0]+=Kpn5*(numa[0][ton1][0]+kua0-numa[0][ton1][1]-kua1)
                        F2[1]-=Kpn5*(numa[0][ton1][0]+kua0-numa[0][ton1][1]-kua1)    
                for ton1 in range(Nma[1]):
                    if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                        F3[0]-=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                        F3[1]+=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                for ton1 in range(Nma1[0]):
                    if dpn==0:
                        kua0=tui2;kua1=0
                    if dpn==1:
                        kua0=0;kua1=tui2  
                    if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                        F1[0]+=Kpn5*(numa1[0][ton1][0]+kua0-numa1[0][ton1][1]-kua1)
                        F1[1]-=Kpn5*(numa1[0][ton1][0]+kua0-numa1[0][ton1][1]-kua1)    
                for ton1 in range(Nma1[1]):#tui2
                      
                    if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                        F4[0]-=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                        F4[1]+=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])  

                Fcc=0
                for x2 in range(pair):
                    if abs(F1[x2])+abs(F2[x2])+abs(F3[x2])+abs(F4[x2])>0.01:
                        Fcc=1
                        break
                if Fcc==0 and abs(F0)<0.0001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
                    de[1][dpn]+=tui1;de[3][dpn]+=tui2
                    for xxx in range(pair):
                        F1[xxx]=p1[xxx]
                        F2[xxx]=p2[xxx]
                        F3[xxx]=p3[xxx]
                        F4[xxx]=p4[xxx]
                    return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]

                jisuan1=[[] for i in range(pair)]
                jisuanp1=[[] for i in range(pair)]
                jisuanpx1=[[] for i in range(pair)]
                for x2 in range(pair):
                    if x2==dpn:
                        for ton1 in range(monum[x2]):
                            jisuan1[x2].append(jisuan[x2][ton1]+tui1)
                        for ton1 in range(monump[x2]):
                            jisuanp1[x2].append(jisuanp[x2][ton1]+tui2)
                        for ton1 in range(monumpx[x2]):
                            jisuanpx1[x2].append(jisuanpx[x2][ton1]-tui1)
                    if x2!=dpn:
                        for ton1 in range(monum[x2]):
                            jisuan1[x2].append(jisuan[x2][ton1])
                        for ton1 in range(monump[x2]):
                            jisuanp1[x2].append(jisuanp[x2][ton1])
                        for ton1 in range(monumpx[x2]):
                            jisuanpx1[x2].append(jisuanpx[x2][ton1])
                #print(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7,'cc1')
                tuf,de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]=jiaocheck(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7)
                if tuf==1:
                    de[1][dpn]+=tui1
                    de[3][dpn]+=tui2
                    for xxx in range(pair):
                        F1[xxx]=p1[xxx]
                        F2[xxx]=p2[xxx]
                        F3[xxx]=p3[xxx]
                        F4[xxx]=p4[xxx]
                    return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]
                    break
    if tq==1:
        lastui1=-1000
        for xn in range(num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]+1):
            if xn < num1[dpn]:
                tui1=(jisuan[dpn][-num1[dpn]+xn]-xdmotor+0.1)
            if xn >= num1[dpn] and xn <num1[dpn]+nump2[dpn]:
                tui1=(jisuanp[dpn][nump2[dpn]-1-xn+num1[dpn]]+xdmotor-0.1)
            if xn >= num1[dpn]+nump2[dpn] and xn <num1[dpn]+nump2[dpn]+num2[dpn]:
                tui1=(jisuan[dpn][num2[dpn]-1-xn+num1[dpn]+nump2[dpn]]+xdmotor-0.1)
            if xn >= num1[dpn]+nump2[dpn]+num2[dpn] and xn <num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]:
                tui1=(jisuanp[dpn][monump[dpn]-nump1[dpn]+xn-num1[dpn]-nump2[dpn]-num2[dpn]]-xdmotor+0.1)
            if xn==num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]:
                tui1=0
            if abs(tui1-lastui1)<0.0001:
                continue
            lastui1=tui1
            lastui2=-1000
            for xn1 in range(numpx2[dpn]+numpx1[dpn]+1):#MT4
                if xn1 <numpx1[dpn]:
                    tui2=(-jisuanpx[dpn][-numpx1[dpn]+xn1]+xdmotor-0.1)
                if xn1 >= numpx1[dpn] and xn1 <numpx1[dpn]+numpx2[dpn]:
                    tui2=-(jisuanpx[dpn][numpx2[dpn]-1-xn1+numpx1[dpn]]+xdmotor-0.1)
                if xn1==numpx1[dpn]+numpx2[dpn]:
                    tui2=0
                if abs(tui2-lastui2)<0.0001 or abs(tui1)+abs(tui2)==0:
                    continue
                lastui2=tui2#####MT3-move-tui1,MT4-move-tui2
                F0=0
                F7=0
                F5=0
                F6=0
                for x2 in range(pair):
                    if x2!=dpn:
                        F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                        F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
                        F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(monum[x2]):
                            if jisuan[x2][ton1]>xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                            if jisuan[x2][ton1]<-xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                        for ton1 in range(monump[x2]):
                            if jisuanp[x2][ton1]>xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                            if jisuanp[x2][ton1]<-xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                        for ton1 in range(monumpx[x2]):
                            if jisuanpx[x2][ton1]>xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                            if jisuanpx[x2][ton1]<-xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                    if x2==dpn:
                        F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                        F7+=Kp1*(l4r[x2]+tui2-Xpo)+Kp2*(l3r[x2]+tui1-Xpo)
                        F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]+tui2-Xkin1)
                        F3[x2]=Kp2*(Xpo-l3r[x2]-tui1)
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2]-tui2)+Kp3*(Xkin1-l4l[x2]-tui2)
                        for ton1 in range(monum[x2]):
                            if jisuan[x2][ton1]-tui1>xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]-tui1-xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]-tui1-xdmotor)*Km
                            if jisuan[x2][ton1]-tui1<-xdmotor:
                                F3[x2]+=(jisuan[x2][ton1]-tui1+xdmotor)*Km
                                F2[x2]+=(jisuan[x2][ton1]-tui1+xdmotor)*Km
                        for ton1 in range(monump[x2]):
                            if jisuanp[x2][ton1]-tui1>xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]-tui1-xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]-tui1-xdmotor)*Km
                            if jisuanp[x2][ton1]-tui1<-xdmotor:
                                F3[x2]+=(jisuanp[x2][ton1]-tui1+xdmotor)*Km
                                F1[x2]+=(jisuanp[x2][ton1]-tui1+xdmotor)*Km
                        for ton1 in range(monumpx[x2]):
                            if jisuanpx[x2][ton1]+tui2>xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]+tui2-xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]+tui2-xdmotor)*Km
                            if jisuanpx[x2][ton1]+tui2<-xdmotor:
                                F4[x2]-=(jisuanpx[x2][ton1]+tui2+xdmotor)*Km
                                F2[x2]-=(jisuanpx[x2][ton1]+tui2+xdmotor)*Km
                for ton1 in range(Nma[0]):
                    if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                        F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                        F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])    
                for ton1 in range(Nma[1]):
                    if dpn==0:
                        kua0=tui1;kua1=0
                    if dpn==1:
                        kua0=0;kua1=tui1    
                    if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                        F3[0]-=Kpn5*(numa[1][ton1][0]+kua0-numa[1][ton1][1]-kua1)
                        F3[1]+=Kpn5*(numa[1][ton1][0]+kua0-numa[1][ton1][1]-kua1)   
                for ton1 in range(Nma1[0]):
                    if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                        F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                        F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])    
                for ton1 in range(Nma1[1]):#tui2
                    if dpn==0:
                        kua0=tui2;kua1=0
                    if dpn==1:
                        kua0=0;kua1=tui2    
                    if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                        F4[0]-=Kpn5*(numa1[1][ton1][0]+kua0-numa1[1][ton1][1]-kua1)
                        F4[1]+=Kpn5*(numa1[1][ton1][0]+kua0-numa1[1][ton1][1]-kua1)  

                Fcc=0
                for x2 in range(pair):
                    if abs(F1[x2])+abs(F2[x2])+abs(F3[x2])+abs(F4[x2])>0.01:
                        Fcc=1
                        break
                if Fcc==0 and abs(F0)<0.0001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
                    de[2][dpn]+=tui1;de[4][dpn]+=tui2
                    for xxx in range(pair):
                        F1[xxx]=p1[xxx]
                        F2[xxx]=p2[xxx]
                        F3[xxx]=p3[xxx]
                        F4[xxx]=p4[xxx]
                    return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]

                jisuan1=[[] for i in range(pair)]
                jisuanp1=[[] for i in range(pair)]
                jisuanpx1=[[] for i in range(pair)]
                for x2 in range(pair):
                    if x2==dpn:
                        for ton1 in range(monum[x2]):
                            jisuan1[x2].append(jisuan[x2][ton1]-tui1)
                        for ton1 in range(monump[x2]):
                            jisuanp1[x2].append(jisuanp[x2][ton1]-tui1)
                        for ton1 in range(monumpx[x2]):
                            jisuanpx1[x2].append(jisuanpx[x2][ton1]+tui2)
                    if x2!=dpn:
                        for ton1 in range(monum[x2]):
                            jisuan1[x2].append(jisuan[x2][ton1])
                        for ton1 in range(monump[x2]):
                            jisuanp1[x2].append(jisuanp[x2][ton1])
                        for ton1 in range(monumpx[x2]):
                            jisuanpx1[x2].append(jisuanpx[x2][ton1])
                #print(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7,'cc2')
                tuf,de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]=jiaocheck(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7)
                if tuf==1:
                    de[2][dpn]+=tui1
                    de[4][dpn]+=tui2
                    for xxx in range(pair):
                        F1[xxx]=p1[xxx]
                        F2[xxx]=p2[xxx]
                        F3[xxx]=p3[xxx]
                        F4[xxx]=p4[xxx]
                    return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]
                    break
    if dpn==0:
        kpn=1
    if dpn==1:
        kpn=0
    lastui1=-1000
    for xn in range(num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]+1):
        if xn < num1[dpn]:
            tui1=(jisuan[dpn][-num1[dpn]+xn]-xdmotor+0.1)
        if xn >= num1[dpn] and xn <num1[dpn]+nump2[dpn]:
            tui1=(jisuanp[dpn][nump2[dpn]-1-xn+num1[dpn]]+xdmotor-0.1)
        if xn >= num1[dpn]+nump2[dpn] and xn <num1[dpn]+nump2[dpn]+num2[dpn]:
            tui1=(jisuan[dpn][num2[dpn]-1-xn+num1[dpn]+nump2[dpn]]+xdmotor-0.1)
        if xn >= num1[dpn]+nump2[dpn]+num2[dpn] and xn <num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]:
            tui1=(jisuanp[dpn][monump[dpn]-nump1[dpn]+xn-num1[dpn]-nump2[dpn]-num2[dpn]]-xdmotor+0.1)
        if xn==num1[dpn]+nump2[dpn]+num2[dpn]+nump1[dpn]:
            tui1=0
        if abs(tui1-lastui1)<0.0001:
            continue
        lastui1=tui1
        lastui2=-1000
        for xn1 in range(numpx2[dpn]+numpx1[dpn]+1):#MT4
            if xn1 <numpx1[dpn]:
                tui2=(-jisuanpx[dpn][-numpx1[dpn]+xn1]+xdmotor-0.1)
            if xn1 >= numpx1[dpn] and xn1 <numpx1[dpn]+numpx2[dpn]:
                tui2=-(jisuanpx[dpn][numpx2[dpn]-1-xn1+numpx1[dpn]]+xdmotor-0.1)
            if xn1==numpx1[dpn]+numpx2[dpn]:
                tui2=0
            if abs(tui2-lastui2)<0.0001:
                continue
            lastui2=tui2
            lastui3=-1000
            for xn2 in range(num1[kpn]+nump2[kpn]+num2[kpn]+nump1[kpn]+1):
                if xn2 < num1[kpn]:
                    tui3=(jisuan[kpn][-num1[kpn]+xn2]-xdmotor+0.1)
                if xn2 >= num1[kpn] and xn2 <num1[kpn]+nump2[kpn]:
                    tui3=(jisuanp[kpn][nump2[kpn]-1-xn2+num1[kpn]]+xdmotor-0.1)
                if xn2 >= num1[kpn]+nump2[kpn] and xn2 <num1[kpn]+nump2[kpn]+num2[kpn]:
                    tui3=(jisuan[kpn][num2[kpn]-1-xn2+num1[kpn]+nump2[kpn]]+xdmotor-0.1)
                if xn2 >= num1[kpn]+nump2[kpn]+num2[kpn] and xn2 <num1[kpn]+nump2[kpn]+num2[kpn]+nump1[kpn]:
                    tui3=(jisuanp[kpn][monump[kpn]-nump1[kpn]+xn2-num1[kpn]-nump2[kpn]-num2[kpn]]-xdmotor+0.1)
                if xn2==num1[kpn]+nump2[kpn]+num2[kpn]+nump1[kpn]:
                    tui3=0
                if abs(tui3-lastui3)<0.0001:
                    continue
                lastui3=tui3
                lastui4=-1000
                for xn3 in range(numpx2[kpn]+numpx1[kpn]+1):#MT4
                    if xn3 <numpx1[kpn]:
                        tui4=(-jisuanpx[kpn][-numpx1[kpn]+xn3]+xdmotor-0.1)
                    if xn3 >= numpx1[kpn] and xn3 <numpx1[kpn]+numpx2[kpn]:
                        tui4=-(jisuanpx[kpn][numpx2[kpn]-1-xn3+numpx1[kpn]]+xdmotor-0.1)
                    if xn3==numpx1[kpn]+numpx2[kpn]:
                        tui4=0
                    if abs(tui4-lastui4)<0.0001 or abs(tui3)+abs(tui4)==0:
                        continue
                    lastui4=tui4
                    #print(tui1,tui2,tui3,tui4)
                    F0=0
                    F7=0
                    F5=0
                    F6=0
                    for x2 in range(pair):
                        if x2!=dpn and x2!=kpn:
                            F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                            F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
                            F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                            F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
                            F3[x2]=Kp2*(Xpo-l3r[x2])
                            F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                            F2[x2]=Kp2*(l2l[x2]-Xpole)
                            F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                            for ton1 in range(monum[x2]):
                                if jisuan[x2][ton1]>xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]-xdmotor)*Km
                                if jisuan[x2][ton1]<-xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]+xdmotor)*Km
                            for ton1 in range(monump[x2]):
                                if jisuanp[x2][ton1]>xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]-xdmotor)*Km
                                if jisuanp[x2][ton1]<-xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]+xdmotor)*Km
                            for ton1 in range(monumpx[x2]):
                                if jisuanpx[x2][ton1]>xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]-xdmotor)*Km
                                if jisuanpx[x2][ton1]<-xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]+xdmotor)*Km
                        if x2==dpn:
                            F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                            F7+=Kp1*(l4r[x2]+tui2-Xpo)+Kp2*(l3r[x2]+tui1-Xpo)
                            F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                            F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]+tui2-Xkin1)
                            F3[x2]=Kp2*(Xpo-l3r[x2]-tui1)
                            F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                            F2[x2]=Kp2*(l2l[x2]-Xpole)
                            F4[x2]=Kp1*(Xpo-l4r[x2]-tui2)+Kp3*(Xkin1-l4l[x2]-tui2)
                            for ton1 in range(monum[x2]):
                                if jisuan[x2][ton1]-tui1>xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]-tui1-xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]-tui1-xdmotor)*Km
                                if jisuan[x2][ton1]-tui1<-xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]-tui1+xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]-tui1+xdmotor)*Km
                            for ton1 in range(monump[x2]):
                                if jisuanp[x2][ton1]-tui1>xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]-tui1-xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]-tui1-xdmotor)*Km
                                if jisuanp[x2][ton1]-tui1<-xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]-tui1+xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]-tui1+xdmotor)*Km
                            for ton1 in range(monumpx[x2]):
                                if jisuanpx[x2][ton1]+tui2>xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]+tui2-xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]+tui2-xdmotor)*Km
                                if jisuanpx[x2][ton1]+tui2<-xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]+tui2+xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]+tui2+xdmotor)*Km
                        if x2==kpn:
                            F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
                            F7+=Kp1*(l4r[x2]+tui4-Xpo)+Kp2*(l3r[x2]+tui3-Xpo)
                            F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
                            F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]+tui4-Xkin1)
                            F3[x2]=Kp2*(Xpo-l3r[x2]-tui3)
                            F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                            F2[x2]=Kp2*(l2l[x2]-Xpole)
                            F4[x2]=Kp1*(Xpo-l4r[x2]-tui4)+Kp3*(Xkin1-l4l[x2]-tui4)
                            for ton1 in range(monum[x2]):
                                if jisuan[x2][ton1]-tui3>xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]-tui3-xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]-tui3-xdmotor)*Km
                                if jisuan[x2][ton1]-tui3<-xdmotor:
                                    F3[x2]+=(jisuan[x2][ton1]-tui3+xdmotor)*Km
                                    F2[x2]+=(jisuan[x2][ton1]-tui3+xdmotor)*Km
                            for ton1 in range(monump[x2]):
                                if jisuanp[x2][ton1]-tui3>xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]-tui3-xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]-tui3-xdmotor)*Km
                                if jisuanp[x2][ton1]-tui3<-xdmotor:
                                    F3[x2]+=(jisuanp[x2][ton1]-tui3+xdmotor)*Km
                                    F1[x2]+=(jisuanp[x2][ton1]-tui3+xdmotor)*Km
                            for ton1 in range(monumpx[x2]):
                                if jisuanpx[x2][ton1]+tui4>xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]+tui4-xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]+tui4-xdmotor)*Km
                                if jisuanpx[x2][ton1]+tui4<-xdmotor:
                                    F4[x2]-=(jisuanpx[x2][ton1]+tui4+xdmotor)*Km
                                    F2[x2]-=(jisuanpx[x2][ton1]+tui4+xdmotor)*Km
                    for ton1 in range(Nma[0]):
                        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                            F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                            F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                    for ton1 in range(Nma[1]):
                        if dpn==0:
                            kua0=tui1;kua1=tui3
                        if dpn==1:
                            kua0=tui3;kua1=tui1    
                        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                            F3[0]-=Kpn5*(numa[1][ton1][0]+kua0-numa[1][ton1][1]-kua1)
                            F3[1]+=Kpn5*(numa[1][ton1][0]+kua0-numa[1][ton1][1]-kua1)  
                    for ton1 in range(Nma1[0]):
                        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                            F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                            F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                    for ton1 in range(Nma1[1]):
                        if dpn==0:
                            kua0=tui2;kua1=tui4
                        if dpn==1:
                            kua0=tui4;kua1=tui2    
                        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                            F4[0]-=Kpn5*(numa1[1][ton1][0]+kua0-numa1[1][ton1][1]-kua1)
                            F4[1]+=Kpn5*(numa1[1][ton1][0]+kua0-numa1[1][ton1][1]-kua1) 
                    Fcc=0
                    for x2 in range(pair):
                        if abs(F1[x2])+abs(F2[x2])+abs(F3[x2])+abs(F4[x2])>0.01:
                            Fcc=1
                            break
                    if Fcc==0 and abs(F0)<0.0001 and abs(F5)<0.00001 and abs(F6)<0.00001 and abs(F7)<0.00001:
                        de[2][dpn]+=tui1;de[4][dpn]+=tui2
                        de[2][kpn]+=tui3;de[4][kpn]+=tui4
                        for xxx in range(pair):
                            F1[xxx]=p1[xxx]
                            F2[xxx]=p2[xxx]
                            F3[xxx]=p3[xxx]
                            F4[xxx]=p4[xxx]
                        return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]

                    jisuan1=[[] for i in range(pair)]
                    jisuanp1=[[] for i in range(pair)]
                    jisuanpx1=[[] for i in range(pair)]
                    for x2 in range(pair):
                        if x2==dpn:
                            for ton1 in range(monum[x2]):
                                jisuan1[x2].append(jisuan[x2][ton1]-tui1)
                            for ton1 in range(monump[x2]):
                                jisuanp1[x2].append(jisuanp[x2][ton1]-tui1)
                            for ton1 in range(monumpx[x2]):
                                jisuanpx1[x2].append(jisuanpx[x2][ton1]+tui2)
                        if x2==kpn:
                            for ton1 in range(monum[x2]):
                                jisuan1[x2].append(jisuan[x2][ton1]-tui3)
                            for ton1 in range(monump[x2]):
                                jisuanp1[x2].append(jisuanp[x2][ton1]-tui3)
                            for ton1 in range(monumpx[x2]):
                                jisuanpx1[x2].append(jisuanpx[x2][ton1]+tui4)
                        if x2!=dpn and x2!=kpn:
                            for ton1 in range(monum[x2]):
                                jisuan1[x2].append(jisuan[x2][ton1])
                            for ton1 in range(monump[x2]):
                                jisuanp1[x2].append(jisuanp[x2][ton1])
                            for ton1 in range(monumpx[x2]):
                                jisuanpx1[x2].append(jisuanpx[x2][ton1])
                    #print(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7,'dd')
                    tuf,de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]=jiaocheck(jisuan1,jisuanp1,jisuanpx1,F0,F1,F2,F3,F4,F5,F6,F7)
                    if tuf==1:
                        de[2][dpn]+=tui1
                        de[4][dpn]+=tui2
                        de[2][kpn]+=tui3;de[4][kpn]+=tui4
                        for xxx in range(pair):
                            F1[xxx]=p1[xxx]
                            F2[xxx]=p2[xxx]
                            F3[xxx]=p3[xxx]
                            F4[xxx]=p4[xxx]
                        return de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]


# In[5]:


def gaussian_elimination(A, b):
    #print(A,b)
    n = len(A)
 
    # Convert the data type of the array to float64.

    A = A.astype(np.float64)
    b = b.astype(np.float64)
    #print(A,b)
    # Gaussian elimination
    for i in range(n - 1):
        max_idx = i
 
        # Pivoting on the column
        for j in range(i + 1, n):
            if abs(A[j][i]) > abs(A[max_idx][i]):
                max_idx = j
 
        # exchange the lines
        A[[i, max_idx]] = A[[max_idx, i]]
        b[[i, max_idx]] = b[[max_idx, i]]
 
        for j in range(i + 1, n):

            multiplier = A[j][i] / A[i][i]
 
            #Update the matrix
            A[j][i:] -= multiplier * A[i][i:]
            b[j] -= multiplier * b[i]
 
    #  solute

    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if A[i][i]!=0 and abs(A[i][i])>1e-12:
            x[i] = (b[i] - np.dot(A[i][i + 1:], x[i + 1:])) / A[i][i]
 
    return x


# In[6]:


def jiaocheck(jisuan,jisuanp,jisuanpx,F0,F1,F2,F3,F4,F5,F6,F7):#Auxiliary Calculation of Displacement
    num3=[0 for i in range(pair)]
    num4=[0 for i in range(pair)]
    nump3=[0 for i in range(pair)]
    nump4=[0 for i in range(pair)]
    numpx3=[0 for i in range(pair)]
    numpx4=[0 for i in range(pair)]
    ji=[[] for i in range(pair)]
    fen0=[[] for i in range(pair)]
    fen1=[[] for i in range(pair)]
    fen2=[[] for i in range(pair)]
    fenp0=[[] for i in range(pair)]
    fenp1=[[] for i in range(pair)]
    fenp2=[[] for i in range(pair)]
    fenpx0=[[] for i in range(pair)]
    fenpx1=[[] for i in range(pair)]
    fenpx2=[[] for i in range(pair)]
    for x2 in range(pair):
        for ton1 in range(monum[x2]):
            if jisuan[x2][ton1]>xdmotor:
                num3[x2]+=1
            if jisuan[x2][ton1]<-xdmotor:
                num4[x2]+=1
        for ton1 in range(monump[x2]):
            if jisuanp[x2][ton1]>xdmotor:
                nump3[x2]+=1
            if jisuanp[x2][ton1]<-xdmotor:
                nump4[x2]+=1
        for ton1 in range(monumpx[x2]):
            if jisuanpx[x2][ton1]>xdmotor:
                numpx3[x2]+=1
            if jisuanpx[x2][ton1]<-xdmotor:
                numpx4[x2]+=1
        m1=[monum[x2]-num3[x2]]
        m0=[num4[x2]]
        m2=[0]
        if monum[x2]>0:
            ttt=jisuan[x2][0]
            for tt1 in range(monum[x2]):
                if abs(jisuan[x2][tt1]-ttt)>0.0001:
                    if jisuan[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuan[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuan[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuan[x2][tt1]
            m2.append(num4[x2])
            m0.append(monum[x2]-num3[x2])
            m1.append(monum[x2])
        fen2[x2]=list(set(m2))
        fen1[x2]=list(set(m1))
        fen0[x2]=list(set(m0))
        fen2[x2].sort(key=m2.index)
        fen1[x2].sort(key=m1.index)
        fen0[x2].sort(key=m0.index)

        m1=[monump[x2]-nump3[x2]]
        m0=[nump4[x2]]
        m2=[0]
        if monump[x2]>0:
            ttt=jisuanp[x2][0]
            for tt1 in range(monump[x2]):
                if abs(jisuanp[x2][tt1]-ttt)>0.0001:
                    if jisuanp[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuanp[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuanp[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuanp[x2][tt1]
            m2.append(nump4[x2])
            m0.append(monump[x2]-nump3[x2])
            m1.append(monump[x2])
        fenp2[x2]=list(set(m2))
        fenp1[x2]=list(set(m1))
        fenp0[x2]=list(set(m0))
        fenp2[x2].sort(key=m2.index)
        fenp1[x2].sort(key=m1.index)
        fenp0[x2].sort(key=m0.index)

        m1=[monumpx[x2]-numpx3[x2]]
        m0=[numpx4[x2]]
        m2=[0]
        if monumpx[x2]>0:
            ttt=jisuanpx[x2][0]
            for tt1 in range(monumpx[x2]):
                if abs(jisuanpx[x2][tt1]-ttt)>0.0001:
                    if jisuanpx[x2][tt1]<-xdmotor:
                        m2.append(tt1)
                    if jisuanpx[x2][tt1]>=-xdmotor and ttt<xdmotor:
                        m0.append(tt1)
                    if jisuanpx[x2][tt1]>xdmotor:
                        m1.append(tt1)
                    ttt=jisuanpx[x2][tt1]
            m2.append(numpx4[x2])
            m0.append(monumpx[x2]-numpx3[x2])
            m1.append(monumpx[x2])
        fenpx2[x2]=list(set(m2))
        fenpx1[x2]=list(set(m1))
        fenpx0[x2]=list(set(m0))
        fenpx2[x2].sort(key=m2.index)
        fenpx1[x2].sort(key=m1.index)
        fenpx0[x2].sort(key=m0.index)
    #########################################to calculate dzk,dz0
    dzk=[[[0 for xxn2 in range(len(fen0[0]))] for xxn1 in range(2)]]
    dz0=[[[0 for xxn2 in range(len(fen2[0]))],[0 for xxn1 in range(len(fen1[0]))]]]
    dzkp=[[[0 for xxn2 in range(len(fenp0[0]))] for xxn1 in range(2)]]
    dz0p=[[[0 for xxn2 in range(len(fenp1[0]))],[0 for xxn1 in range(len(fenp2[0]))]]]
    dzkpx=[[[0 for xxn2 in range(len(fenpx0[0]))] for xxn1 in range(2)]]
    dz0px=[[[0 for xxn2 in range(len(fenpx1[0]))],[0 for xxn1 in range(len(fenpx2[0]))]]]
    for x3 in range(1,pair):
        dzk.append([[0 for xxn2 in range(len(fen0[x3]))] for xxn1 in range(2)])
        dz0.append([[0 for xxn2 in range(len(fen2[x3]))],[0 for xxn1 in range(len(fen1[x3]))]])
        dzkp.append([[0 for xxn2 in range(len(fenp0[x3]))] for xxn1 in range(2)])
        dz0p.append([[0 for xxn2 in range(len(fenp1[x3]))],[0 for xxn1 in range(len(fenp2[x3]))]])
        dzkpx.append([[0 for xxn2 in range(len(fenpx0[x3]))] for xxn1 in range(2)])
        dz0px.append([[0 for xxn2 in range(len(fenpx1[x3]))],[0 for xxn1 in range(len(fenpx2[x3]))]])
    for x2 in range(pair):
        ta1=-1
        for xck in fen0[x2]:
            ta1+=1
            cnnk=monum[x2]-num3[x2]-fen0[x2][-ta1-1]
            chek=monum[x2]-num3[x2]-1#dxc增，ds1>0
            for xxxn in range(cnnk):
                dzk[x2][0][ta1]+=jisuan[x2][chek-xxxn]-xdmotor
        ta1=-1
        for xck in fen0[x2]:
            ta1+=1
            cnnk=xck-num4[x2]
            chek=num4[x2]#dxc减，ds1<0
            for xxxn in range(cnnk):
                dzk[x2][1][ta1]+=jisuan[x2][chek+xxxn]+xdmotor 
        ta2=-1
        for xckp in fenp0[x2]:
            ta2+=1
            cnnkp=xckp-nump4[x2]
            chekp=nump4[x2]#0对应dxcp减
            for xxxn in range(cnnkp):
                dzkp[x2][0][ta2]+=jisuanp[x2][chekp+xxxn]+xdmotor
        ta2=-1
        for xckp in fenp0[x2]:
            ta2+=1
            cnnkp=monump[x2]-nump3[x2]-fenp0[x2][-ta2-1]
            chekp=monump[x2]-nump3[x2]-1 
            for xxxn in range(cnnkp):
                dzkp[x2][1][ta2]+=jisuanp[x2][chekp-xxxn]-xdmotor 
        ta2=-1
        for xckpx in fenpx0[x2]:
            ta2+=1
            cnnkpx=xckpx-numpx4[x2]
            chekpx=numpx4[x2]
            for xxxn in range(cnnkpx):
                dzkpx[x2][0][ta2]+=jisuanpx[x2][chekpx+xxxn]+xdmotor
        ta2=-1
        for xckpx in fenpx0[x2]:
            ta2+=1
            cnnkpx=monumpx[x2]-numpx3[x2]-fenpx0[x2][-ta2-1]
            chekpx=monumpx[x2]-numpx3[x2]-1 
            for xxxn in range(cnnkpx):
                dzkpx[x2][1][ta2]+=jisuanpx[x2][chekpx-xxxn]-xdmotor 
        ta3=-1
        for xc0 in fen2[x2]:
            ta3+=1
            cnn0=num4[x2]-fen2[x2][-ta3-1]
            che0=num4[x2]-1
            for xxxn in range(cnn0):
                dz0[x2][0][ta3]+=jisuan[x2][che0-xxxn]+xdmotor
        ta3=-1
        for xc0 in fen1[x2]:    
            ta3+=1
            cnn0=xc0-monum[x2]+num3[x2]
            che0=monum[x2]-num3[x2]
            for xxxn in range(cnn0):
                dz0[x2][1][ta3]+=jisuan[x2][che0+xxxn]-xdmotor
        ta4=-1
        for xc0p in fenp1[x2]:
            ta4+=1
            cnn0p=xc0p-monump[x2]+nump3[x2]
            che0p=monump[x2]-nump3[x2]
            for xxxn in range(cnn0p):
                dz0p[x2][0][ta4]+=jisuanp[x2][che0p+xxxn]-xdmotor
        ta4=-1
        for xc0p in fenp2[x2]:
            ta4+=1
            cnn0p=nump4[x2]-fenp2[x2][-ta4-1]
            che0p=nump4[x2]-1
            for xxxn in range(cnn0p):
                dz0p[x2][1][ta4]+=jisuanp[x2][che0p-xxxn]+xdmotor
        ta4=-1
        for xc0px in fenpx1[x2]:
            ta4+=1
            cnn0px=xc0px-monumpx[x2]+numpx3[x2]
            che0px=monumpx[x2]-numpx3[x2]
            for xxxn in range(cnn0px):
                dz0px[x2][0][ta4]+=jisuanpx[x2][che0px+xxxn]-xdmotor
        ta4=-1
        for xc0px in fenpx2[x2]:
            ta4+=1
            cnn0px=numpx4[x2]-fenpx2[x2][-ta4-1]
            che0px=numpx4[x2]-1
            for xxxn in range(cnn0px):
                dz0px[x2][1][ta4]+=jisuanpx[x2][che0px-xxxn]+xdmotor
   
    #print(round(F0),round(F1,2),round(F2,2),round(F3,2),round(F4,2),round(F5,2),round(F6,2),round(F7,2))
    for x2 in range(pair):
        for xx1 in range(2):
            ta1=-1
            for xck in fen0[x2]:
                ta1+=1
                if xx1==0: 
                    if ta1<len(dzk[x2][0]):
                        cnnk=monum[x2]-num3[x2]-fen0[x2][-ta1-1]
                    if ta1>=len(dzk[x2][0]):
                        break
                if xx1==1:
                    if ta1<len(dzk[x2][1]):
                        cnnk=xck-num4[x2]
                    if ta1>=len(dzk[x2][1]):
                        break
                for xx2 in range(2):
                    ta2=-1
                    for xckp in fenp0[x2]:
                        ta2+=1
                        if xx2==0:
                            if ta2<len(dzkp[x2][0]):
                                cnnkp=xckp-nump4[x2]
                            if ta2>=len(dzkp[x2][0]):
                                break
                        if xx2==1:
                            if ta2<len(dzkp[x2][1]):
                                cnnkp=monump[x2]-nump3[x2]-fenp0[x2][-ta2-1]
                            if ta2>=len(dzkp[x2][1]):
                                break
                        for xx3 in range(2):
                            ta5=-1
                            for xckpx in fenpx0[x2]:    
                                ta5+=1
                                if xx3==0:
                                    if ta5<len(dzkpx[x2][0]):
                                        cnnkpx=xckpx-numpx4[x2]
                                    if ta5>=len(dzkpx[x2][0]):
                                        break
                                if xx3==1:
                                    if ta5<len(dzkpx[x2][1]):
                                        cnnkpx=monumpx[x2]-numpx3[x2]-fenpx0[x2][-ta5-1]
                                    if ta5>=len(dzkpx[x2][1]):
                                        break
                                if xx1==0:
                                    xa3=fen2[x2]
                                if xx1==1:
                                    xa3=fen1[x2]
                                ta3=-1
                                for xc0 in xa3:
                                    #print(xqq,xc0,len(dz0[xqq][1]),len(dz0[xqq][0]),'3')
                                    ta3+=1  
                                    if xx1==0:
                                        if ta3<len(dz0[x2][0]):
                                            cnn0=num4[x2]-fen2[x2][-ta3-1]
                                        if ta3>=len(dz0[x2][0]):
                                            break
                                    if xx1==1:
                                        if ta3<len(dz0[x2][1]):
                                            cnn0=xc0-monum[x2]+num3[x2]
                                        if ta3>=len(dz0[x2][1]):
                                            break
                                    if xx2==0:
                                        xa4=fenp1[x2]
                                    if xx2==1:
                                        xa4=fenp2[x2]
                                    ta4=-1
                                    for xc0p in xa4:
                                        ta4+=1
                                        if xx2==0:
                                            if ta4<len(dz0p[x2][0]):
                                                cnn0p=xc0p-monump[x2]+nump3[x2]
                                            if ta4>=len(dz0p[x2][0]):
                                                break
                                        if xx2==1:
                                            if ta4<len(dz0p[x2][1]):
                                                cnn0p=nump4[x2]-fenp2[x2][-ta4-1]    
                                        if xx3==0:
                                            xa5=fenpx1[x2]
                                        if xx3==1:
                                            xa5=fenpx2[x2]
                                        ta6=-1
                                        for xc0px in xa5:
                                            ta6+=1
                                            if xx3==0:
                                                if ta6<len(dz0px[x2][0]):
                                                    cnn0px=xc0px-monumpx[x2]+numpx3[x2]
                                                if ta6>=len(dz0px[x2][0]):
                                                    break
                                            if xx3==1:
                                                if ta6<len(dz0px[x2][1]):
                                                    cnn0px=numpx4[x2]-fenpx2[x2][-ta6-1] 
                                                if ta6>len(dz0px[x2][1]):
                                                    break
                                            if cnn0+cnnk==0 and xx1!=0:
                                                continue
                                            if cnn0p+cnnkp==0 and xx2!=0:
                                                continue
                                            if cnn0px+cnnkpx==0 and xx3!=0:
                                                continue    
                                            
                                            ji[x2].append([cnn0,cnnk,cnn0p,cnnkp,cnn0px,cnnkpx,xx1,xx2,xx3,ta1,ta2,ta3,ta4,ta5,ta6])
    a00=round(2*(-Kp1-Kp2),3);a10=Kp2;a30=Kp1;a80=Kp2;a100=Kp1
    a21=Kp2;a41=Kp1;a71=round((-Kp1-Kp2)*2,3);a91=Kp2;a111a=Kp1
    a32=Kp3;a52=round((-Kp3-Kp4)*2,3);a62=round(Kp4*2,3);a102=Kp3
    a43=Kp3;a53=round(Kp4*2,3);a63=round((-Kp3-Kp4 )*2,3);a113=Kp3
    zq=[0 for i in range(pair)]
    cxx=[0 for i in range(pair)]
    for cx0 in range(len(ji[0])):
        cxx[0]=cx0
        if abs(sum(zq)-pair)<0.01:
            break
        for cx1 in range(len(ji[1])):
            cxx[1]=cx1
            if abs(sum(zq)-pair)<0.01:
                break

            cnn0=ji[0][cxx[0]][0]
            cnnk=ji[0][cxx[0]][1]
            cnn0p=ji[0][cxx[0]][2]
            cnnkp=ji[0][cxx[0]][3]
            cnn0px=ji[0][cxx[0]][4]
            cnnkpx=ji[0][cxx[0]][5]
            xx1=ji[0][cxx[0]][6]
            xx2=ji[0][cxx[0]][7]
            xx3=ji[0][cxx[0]][8]
            ta1=ji[0][cxx[0]][9]
            ta2=ji[0][cxx[0]][10]
            ta3=ji[0][cxx[0]][11]
            ta4=ji[0][cxx[0]][12]
            ta5=ji[0][cxx[0]][13]
            ta6=ji[0][cxx[0]][14]
            a04=-Kp1;a24=round(-(nump3[0]+nump4[0]+cnnkp-cnn0p)*Km,3);a34=round((nump3[0]+nump4[0]+cnnkp-cnn0p)*Km+Kp1+Kp3+Kpn5*nunum1[0],3)
            a54=-Kp3;a104=round(-Kpn5*nunum1[0],3)
            a05=-Kp2;a15=round(((numpx3[0]+numpx4[0]+cnnkpx-cnn0px)+(num3[0]+num4[0]+cnnk-cnn0))*Km+Kp2+Kpn5*nunum[0],3);a85=round(-Kpn5*nunum[0],3)
            a25=round(-Km*(num3[0]+num4[0]+cnnk-cnn0),3);a45=round(-(numpx3[0]+numpx4[0]+cnnkpx-cnn0px)*Km,3)
            a16=round((num3[0]+num4[0]+cnnk-cnn0)*Km,3);a26=round(-Km*(nump3[0]+nump4[0]+cnnkp-cnn0p+num3[0]+num4[0]+cnnk-cnn0)-Kp2-Kpn5*nunum[1],3)
            a36=round(Km*(nump3[0]+nump4[0]+cnnkp-cnn0p),3);a76=Kp2;a96=round(Kpn5*nunum[1],3)
            a17=round(Km*(numpx3[0]+numpx4[0]+cnnkpx-cnn0px),3);a47=round(-Km*(numpx3[0]+numpx4[0]+cnnkpx-cnn0px)-Kp1-Kp3-Kpn5*nunum1[1],3)
            a67=Kp3;a77=Kp1;a117=round(Kpn5*nunum1[1],3)
            b4=-Km*(dzkp[0][xx2][ta2]-dz0p[0][xx2][ta4])-F1[0]
            b5=-Km*(-dzkpx[0][xx3][ta5]+dz0px[0][xx3][ta6]+dzk[0][xx1][ta1]-dz0[0][xx1][ta3])-F2[0]
            b6=-Km*(dzkp[0][xx2][ta2]-dz0p[0][xx2][ta4]+dzk[0][xx1][ta1]-dz0[0][xx1][ta3])-F3[0]
            b7=Km*(dzkpx[0][xx3][ta5]-dz0px[0][xx3][ta6])-F4[0]
            cnn0=ji[1][cxx[1]][0]
            cnnk=ji[1][cxx[1]][1]
            cnn0p=ji[1][cxx[1]][2]
            cnnkp=ji[1][cxx[1]][3]
            cnn0px=ji[1][cxx[1]][4]
            cnnkpx=ji[1][cxx[1]][5]
            xx1=ji[1][cxx[1]][6]
            xx2=ji[1][cxx[1]][7]
            xx3=ji[1][cxx[1]][8]
            ta1=ji[1][cxx[1]][9]
            ta2=ji[1][cxx[1]][10]
            ta3=ji[1][cxx[1]][11]
            ta4=ji[1][cxx[1]][12]
            ta5=ji[1][cxx[1]][13]
            ta6=ji[1][cxx[1]][14]
            a08=-Kp1;a28=round(-(nump3[1]+nump4[1]+cnnkp-cnn0p)*Km,3);a38=round((nump3[1]+nump4[1]+cnnkp-cnn0p)*Km+Kp1+Kp3+Kpn5*nunum1[0],3)
            a58=-Kp3;a038=round(-Kpn5*nunum1[0],3)
            a09=-Kp2;a19=round(((numpx3[1]+numpx4[1]+cnnkpx-cnn0px)+(num3[1]+num4[1]+cnnk-cnn0))*Km+Kp2+Kpn5*nunum[0],3);a019=round(-Kpn5*nunum[0],3)
            a29=round(-Km*(num3[1]+num4[1]+cnnk-cnn0),3);a49=round(-(numpx3[1]+numpx4[1]+cnnkpx-cnn0px)*Km,3)
            a110=round((num3[1]+num4[1]+cnnk-cnn0)*Km,3);a210=round(-Km*(nump3[1]+nump4[1]+cnnkp-cnn0p+num3[1]+num4[1]+cnnk-cnn0)-Kp2-Kpn5*nunum[1],3)
            a310=round(Km*(nump3[1]+nump4[1]+cnnkp-cnn0p),3);a710=Kp2;a0210=round(Kpn5*nunum[1],3)
            a111=round(Km*(numpx3[1]+numpx4[1]+cnnkpx-cnn0px),3);a411=round(-Km*(numpx3[1]+numpx4[1]+cnnkpx-cnn0px)-Kp1-Kp3-Kpn5*nunum1[1],3)
            a611=Kp3;a711=Kp1;a0411=round(Kpn5*nunum1[1],3)
            b8=-Km*(dzkp[1][xx2][ta2]-dz0p[1][xx2][ta4])-F1[1]
            b9=-Km*(-dzkpx[1][xx3][ta5]+dz0px[1][xx3][ta6]+dzk[1][xx1][ta1]-dz0[1][xx1][ta3])-F2[1]
            b10=-Km*(dzkp[1][xx2][ta2]-dz0p[1][xx2][ta4]+dzk[1][xx1][ta1]-dz0[1][xx1][ta3])-F3[1]
            b11=Km*(dzkpx[1][xx3][ta5]-dz0px[1][xx3][ta6])-F4[1] 

            b=np.array([round(-F0,14),round(-F7,14),round(-F5,14),round(-F6,14),round(b4,14),round(b5,14),round(b6,14),round(b7,14),round(b8,14),round(b9,14),round(b10,14),round(b11,14)])
            try:
                A=np.array([[a00,a10,0,a30,0,0,0,0,a80,0,a100,0],[0,0,a21,0,a41,0,0,a71,0,a91,0,a111a],
                                [0,0,0,a32,0,a52,a62,0,0,0,a102,0],[0,0,0,0,a43,a53,a63,0,0,0,0,a113],
                                [a04,0,a24,a34,0,a54,0,0,0,0,a104,0],[a05,a15,a25,0,a45,0,0,0,a85,0,0,0],
                                [0,a16,a26,a36,0,0,0,a76,0,a96,0,0],[0,a17,0,0,a47,0,a67,a77,0,0,0,a117],
                                [a08,0,0,a038,0,a58,0,0,0,a28,a38,0],[a09,a019,0,0,0,0,0,0,a19,a29,0,a49],
                                [0,0,a0210,0,0,0,0,a710,a110,a210,a310,0],[0,0,0,0,a0411,0,a611,a711,a111,0,0,a411]])
                #b=np.array([-F0,-F7,-F5,-F6,b4,b5,b6,b7])

                mtui=np.linalg.solve(A,b)

            except:
                #print(A,b)
                mtui = gaussian_elimination(A,b)
               
            tuzx=(mtui[0]+mtui[7])/2
            for xn in range(12):
                mtui[xn]-=tuzx
            de=[mtui[0],[mtui[1],mtui[8]],[mtui[2],mtui[9]],[mtui[3],mtui[10]],[mtui[4],mtui[11]],mtui[5],mtui[6],mtui[7]]
            #print(de)
            chezq=0
            for x3 in range(pair):
                check0=-1#################from K to 0
                checkk=-1################from 0 to K
                checkp0=-1#################from K to 0
                checkpk=-1#################from 0 to K
                checkpx0=-1###############
                checkpxk=-1################
                cnn0=ji[x3][cxx[x3]][0]
                cnnk=ji[x3][cxx[x3]][1]
                cnn0p=ji[x3][cxx[x3]][2]
                cnnkp=ji[x3][cxx[x3]][3]
                cnn0px=ji[x3][cxx[x3]][4]
                cnnkpx=ji[x3][cxx[x3]][5]
                xx1=ji[x3][cxx[x3]][6]
                xx2=ji[x3][cxx[x3]][7]
                xx3=ji[x3][cxx[x3]][8]
                #print(chezq)
                if cnn0+cnnk==0:
                    for ton1 in range(monum[x3]):
                        if abs(jisuan[x3][ton1])<=xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuan[x3][ton1])>xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0+cnnk!=0:    
                    if xx1==0 and de[1][x3]-de[2][x3]>=0:
                        check0=0###############
                        checkk=0###############
                        for ton1 in range(monum[x3]):
                            if jisuan[x3][ton1]<-xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                                check0+=1
                            if abs(jisuan[x3][ton1])<=xdmotor and jisuan[x3][ton1]+de[1][x3]-de[2][x3]>xdmotor:
                                checkk+=1
                    if xx1==1 and de[1][x3]-de[2][x3]<=0:
                        check0=0################
                        checkk=0################
                        for ton1 in range(monum[x3]):
                            if abs(jisuan[x3][ton1])<=xdmotor and jisuan[x3][ton1]+de[1][x3]-de[2][x3]<-xdmotor:
                                checkk+=1
                            if jisuan[x3][ton1]>xdmotor and abs(jisuan[x3][ton1]+de[1][x3]-de[2][x3])<=xdmotor:
                                check0+=1
                    if check0!=cnn0 or checkk!=cnnk:
                        chezq=-1
                        break
                #print(chezq)
                if chezq==-1:
                    #print(x3)
                    break
                if cnn0p+cnnkp==0:
                    for ton1 in range(monump[x3]):
                        if abs(jisuanp[x3][ton1])<=xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuanp[x3][ton1])>xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0p+cnnkp!=0:
                    if xx2==1 and de[3][x3]-de[2][x3]>=0:
                        checkp0=0################
                        checkpk=0###############
                        for ton1 in range(monump[x3]):
                            if jisuanp[x3][ton1]<-xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                                checkp0+=1
                            if abs(jisuanp[x3][ton1])<=xdmotor and jisuanp[x3][ton1]+de[3][x3]-de[2][x3]>xdmotor:
                                checkpk+=1
                    if xx2==0 and de[3][x3]-de[2][x3]<=0:#dxcp减
                        checkp0=0
                        checkpk=0
                        for ton1 in range(monump[x3]):
                            if abs(jisuanp[x3][ton1])<=xdmotor and jisuanp[x3][ton1]+de[3][x3]-de[2][x3]<-xdmotor:
                                checkpk+=1
                            if jisuanp[x3][ton1]>xdmotor and abs(jisuanp[x3][ton1]+de[3][x3]-de[2][x3])<=xdmotor:
                                checkp0+=1
                    if checkp0!=cnn0p or checkpk!=cnnkp:
                        chezq=-1
                        break
                #print(chezq)
                if chezq==-1:
                    #print(x3,'p')
                    break
                if cnn0px+cnnkpx==0:
                    for ton1 in range(monumpx[x3]):
                        if abs(jisuanpx[x3][ton1])<=xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])>xdmotor:
                            chezq=-1
                            break
                        if abs(jisuanpx[x3][ton1])>xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                            chezq=-1
                            break
                #print(chezq)
                if cnn0px+cnnkpx!=0:
                    if xx3==1 and de[4][x3]-de[1][x3]>=0:
                        checkpx0=0###############
                        checkpxk=0################
                        for ton1 in range(monumpx[x3]):
                            if jisuanpx[x3][ton1]<-xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                                checkpx0+=1
                            if abs(jisuanpx[x3][ton1])<=xdmotor and jisuanpx[x3][ton1]+de[4][x3]-de[1][x3]>xdmotor:
                                checkpxk+=1
                    if xx3==0 and de[4][x3]-de[1][x3]<=0:#
                        checkpx0=0#################
                        checkpxk=0#################
                        for ton1 in range(monumpx[x3]):
                            if abs(jisuanpx[x3][ton1])<=xdmotor and jisuanpx[x3][ton1]+de[4][x3]-de[1][x3]<-xdmotor:
                                checkpxk+=1
                            if jisuanpx[x3][ton1]>xdmotor and abs(jisuanpx[x3][ton1]+de[4][x3]-de[1][x3])<=xdmotor:
                                checkpx0+=1
                    if checkpx0!=cnn0px or checkpxk!=cnnkpx:
                        chezq=-1
                        break

                if chezq==-1:
                    break
            if chezq==-1:
                continue
            return 1,de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]
    return 0,de[0],de[1],de[2],de[3],de[4],de[5],de[6],de[7]





def energy(x1,i):#calculate the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by two iMTs taking a forward step and a backward step

    ypsl11=0
    ypsl22=0
    for x2 in range(pair):
        ypsl11+=Kp3*(((Xkin+zj5[x1][i]-l1r[x2]-zj3[x1][i][x2])**2+(Xkin1+zj6[x1][i]-l4l[x2]-zj4[x1][i][x2])**2)/2)+Kp2*((Xpole+zj0[x1][i]-l2l[x2]-zj1[x1][i][x2])**2+(Xpo+zj7[x1][i]-l3r[x2]-zj2[x1][i][x2])**2)/2+Kp1*((Xpole+zj0[x1][i]-l1l[x2]-zj3[x1][i][x2])**2+(Xpo+zj7[x1][i]-l4r[x2]-zj4[x1][i][x2])**2)/2+Kp4*((Xkin1+zj6[x1][i]-Xkin-zj5[x1][i]-kdis)**2)/2
        ypsl22+=Kp3*(((Xkin+js5[x1][i]-l1r[x2]-js3[x1][i][x2])**2+(Xkin1+js6[x1][i]-l4l[x2]-js4[x1][i][x2])**2)/2)+Kp2*((Xpole+js0[x1][i]-l2l[x2]-js1[x1][i][x2])**2+(Xpo+js7[x1][i]-l3r[x2]-js2[x1][i][x2])**2)/2+Kp1*((Xpole+js0[x1][i]-l1l[x2]-js3[x1][i][x2])**2+(Xpo+js7[x1][i]-l4r[x2]-js4[x1][i][x2])**2)/2+Kp4*((Xkin1+js6[x1][i]-Xkin-js5[x1][i]-kdis)**2)/2
            #ypsl11+=Kp3*(((Xkin-l1r[x2])**2+(Xkin1-l4l[x2])**2)/2)+Kp2*((Xpole-l2l[x2])**2+(Xpo-l3r[x2])**2)/2+Kp1*((Xpole-l1l[x2])**2+(Xpo-l4r[x2])**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
            #ypsl22+=Kp3*(((Xkin-l1r[x2])**2+(Xkin1-l4l[x2])**2)/2)+Kp2*((Xpole-l2l[x2])**2+(Xpo-l3r[x2])**2)/2+Kp1*((Xpole-l1l[x2])**2+(Xpo-l4r[x2])**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
        
        for xn in range(nu[x2]):
            if kina[x2][xn]==1 and kinb[x2][xn]==1:
                if xn==i and x2==x1:
                    ppx=dx1+d+zj1[x1][i][x2]-zj2[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=dx1-d+js1[x1][i][x2]-js2[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
                if xn!=i or x2!=x1:
                    ppx=dxc[x2][xn]+zj1[x1][i][x2]-zj2[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=dxc[x2][xn]+js1[x1][i][x2]-js2[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1[x2]):
            if kinap[x2][xn]+kinbp[x2][xn]==2:
                ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+zj3[x1][i][x2]-zj2[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+js3[x1][i][x2]-js2[x1][i][x2]
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1x[x2]):
            if kinapx[x2][xn]+kinbpx[x2][xn]==2:
                ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+zj4[x1][i][x2]-zj1[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+js4[x1][i][x2]-js1[x1][i][x2]
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
    for ton1 in range(Nma[0]):   
        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
            ypsl11+=Kpn5*((numa[0][ton1][0]+zj1[x1][i][0]-numa[0][ton1][1]-zj1[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[0][ton1][0]+js1[x1][i][0]-numa[0][ton1][1]-js1[x1][i][1])**2)/2
    for ton1 in range(Nma[1]):   
        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
            ypsl11+=Kpn5*((numa[1][ton1][0]+zj2[x1][i][0]-numa[1][ton1][1]-zj2[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[1][ton1][0]+js2[x1][i][0]-numa[1][ton1][1]-js2[x1][i][1])**2)/2 
    for ton1 in range(Nma1[0]):   
        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
            ypsl11+=Kpn5*((numa1[0][ton1][0]+zj3[x1][i][0]-numa1[0][ton1][1]-zj3[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[0][ton1][0]+js3[x1][i][0]-numa1[0][ton1][1]-js3[x1][i][1])**2)/2
    for ton1 in range(Nma1[1]):   
        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
            ypsl11+=Kpn5*((numa1[1][ton1][0]+zj4[x1][i][0]-numa1[1][ton1][1]-zj4[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[1][ton1][0]+js4[x1][i][0]-numa1[1][ton1][1]-js4[x1][i][1])**2)/2 
    return ypsl11,ypsl22


# In[9]:


def energyp(x1,i):#calculate the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by the iMT connected to the right spindle pole and the kMT connected to the left spindle pole taking a forward step and a backward step
    ypsl11=0
    ypsl22=0
    for x2 in range(pair):
        ypsl11+=Kp3*(((Xkin+zj5p[x1][i]-l1r[x2]-zj3p[x1][i][x2])**2+(Xkin1+zj6p[x1][i]-l4l[x2]-zj4p[x1][i][x2])**2)/2)+Kp2*((Xpole+zj0p[x1][i]-l2l[x2]-zj1p[x1][i][x2])**2+(Xpo+zj7p[x1][i]-l3r[x2]-zj2p[x1][i][x2])**2)/2+Kp1*((Xpole+zj0p[x1][i]-l1l[x2]-zj3p[x1][i][x2])**2+(Xpo+zj7p[x1][i]-l4r[x2]-zj4p[x1][i][x2])**2)/2+Kp4*((Xkin1+zj6p[x1][i]-Xkin-zj5p[x1][i]-kdis)**2)/2
        ypsl22+=Kp3*(((Xkin+js5p[x1][i]-l1r[x2]-js3p[x1][i][x2])**2+(Xkin1+js6p[x1][i]-l4l[x2]-js4p[x1][i][x2])**2)/2)+Kp2*((Xpole+js0p[x1][i]-l2l[x2]-js1p[x1][i][x2])**2+(Xpo+js7p[x1][i]-l3r[x2]-js2p[x1][i][x2])**2)/2+Kp1*((Xpole+js0p[x1][i]-l1l[x2]-js3p[x1][i][x2])**2+(Xpo+js7p[x1][i]-l4r[x2]-js4p[x1][i][x2])**2)/2+Kp4*((Xkin1+js6p[x1][i]-Xkin-js5p[x1][i]-kdis)**2)/2

        for xn in range(nu[x2]):
            if kina[x2][xn]==1 and kinb[x2][xn]==1:
                ppx=dxc[x2][xn]+zj1p[x1][i][x2]-zj2p[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=dxc[x2][xn]+js1p[x1][i][x2]-js2p[x1][i][x2]    
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1[x2]):
            if xn==i and x2==x1:
                if kinap[x2][xn]+kinbp[x2][xn]==2:
                    ppx=dx1+d+zj3p[x1][i][x2]-zj2p[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=dx1-d+js3p[x1][i][x2]-js2p[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
            if xn!=i or x2!=x1:
                if kinap[x2][xn]+kinbp[x2][xn]==2:
                    ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+zj3p[x1][i][x2]-zj2p[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+js3p[x1][i][x2]-js2p[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1x[x2]):
            if kinapx[x2][xn]+kinbpx[x2][xn]==2:
                ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+zj4p[x1][i][x2]-zj1p[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+js4p[x1][i][x2]-js1p[x1][i][x2]
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
    for ton1 in range(Nma[0]):   
        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
            ypsl11+=Kpn5*((numa[0][ton1][0]+zj1p[x1][i][0]-numa[0][ton1][1]-zj1p[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[0][ton1][0]+js1p[x1][i][0]-numa[0][ton1][1]-js1p[x1][i][1])**2)/2
    for ton1 in range(Nma[1]):   
        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
            ypsl11+=Kpn5*((numa[1][ton1][0]+zj2p[x1][i][0]-numa[1][ton1][1]-zj2p[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[1][ton1][0]+js2p[x1][i][0]-numa[1][ton1][1]-js2p[x1][i][1]) **2)/2
    for ton1 in range(Nma1[0]):   
        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
            ypsl11+=Kpn5*((numa1[0][ton1][0]+zj3p[x1][i][0]-numa1[0][ton1][1]-zj3p[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[0][ton1][0]+js3p[x1][i][0]-numa1[0][ton1][1]-js3p[x1][i][1])**2)/2
    for ton1 in range(Nma1[1]):   
        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
            ypsl11+=Kpn5*((numa1[1][ton1][0]+zj4p[x1][i][0]-numa1[1][ton1][1]-zj4p[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[1][ton1][0]+js4p[x1][i][0]-numa1[1][ton1][1]-js4p[x1][i][1])**2)/2 
    return ypsl11,ypsl22


# In[10]:


def energypx(x1,i):#calculate the energy of the spindle systems after the kinesin-5 motors in the overlap region formed by the iMT connected to the left spindle pole and the kMT connected to the right spindle pole taking a forward step and a backward step 
    ypsl11=0
    ypsl22=0
    for x2 in range(pair):
        ypsl11+=Kp3*(((Xkin+zj5px[x1][i]-l1r[x2]-zj3px[x1][i][x2])**2+(Xkin1+zj6px[x1][i]-l4l[x2]-zj4px[x1][i][x2])**2)/2)+Kp2*((Xpole+zj0px[x1][i]-l2l[x2]-zj1px[x1][i][x2])**2+(Xpo+zj7px[x1][i]-l3r[x2]-zj2px[x1][i][x2])**2)/2+Kp1*((Xpole+zj0px[x1][i]-l1l[x2]-zj3px[x1][i][x2])**2+(Xpo+zj7px[x1][i]-l4r[x2]-zj4px[x1][i][x2])**2)/2+Kp4*((Xkin1+zj6px[x1][i]-Xkin-zj5px[x1][i]-kdis)**2)/2
        ypsl22+=Kp3*(((Xkin+js5px[x1][i]-l1r[x2]-js3px[x1][i][x2])**2+(Xkin1+js6px[x1][i]-l4l[x2]-js4px[x1][i][x2])**2)/2)+Kp2*((Xpole+js0px[x1][i]-l2l[x2]-js1px[x1][i][x2])**2+(Xpo+js7px[x1][i]-l3r[x2]-js2px[x1][i][x2])**2)/2+Kp1*((Xpole+js0px[x1][i]-l1l[x2]-js3px[x1][i][x2])**2+(Xpo+js7px[x1][i]-l4r[x2]-js4px[x1][i][x2])**2)/2+Kp4*((Xkin1+js6px[x1][i]-Xkin-js5px[x1][i]-kdis)**2)/2
        for xn in range(nu[x2]):
            if kina[x2][xn]==1 and kinb[x2][xn]==1:
                ppx=dxc[x2][xn]+zj1px[x1][i][x2]-zj2px[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=dxc[x2][xn]+js1px[x1][i][x2]-js2px[x1][i][x2]
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1[x2]):
            if kinap[x2][xn]+kinbp[x2][xn]==2:
                ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+zj3px[x1][i][x2]-zj2px[x1][i][x2]
                if ppx>xdmotor:
                    ypsl11+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl11+=Km*((ppx+xdmotor)**2)/2
                ppx=motorp[x2][xn][0]-motorp[x2][xn][1]+js3px[x1][i][x2]-js2px[x1][i][x2]
                if ppx>xdmotor:
                    ypsl22+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl22+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1x[x2]):
            if xn==i and x2==x1:
                if kinapx[x2][xn]+kinbpx[x2][xn]==2:
                    ppx=dx1+d+zj4px[x1][i][x2]-zj1px[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=dx1-d+js4px[x1][i][x2]-js1px[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
            if xn!=i or x2!=x1:
                if kinapx[x2][xn]+kinbpx[x2][xn]==2:
                    ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+zj4px[x1][i][x2]-zj1px[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl11+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl11+=Km*((ppx+xdmotor)**2)/2
                    ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]+js4px[x1][i][x2]-js1px[x1][i][x2]
                    if ppx>xdmotor:
                        ypsl22+=Km*((ppx-xdmotor)**2)/2
                    if ppx<-xdmotor:
                        ypsl22+=Km*((ppx+xdmotor)**2)/2
    for ton1 in range(Nma[0]):   
        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
            ypsl11+=Kpn5*((numa[0][ton1][0]+zj1px[x1][i][0]-numa[0][ton1][1]-zj1px[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[0][ton1][0]+js1px[x1][i][0]-numa[0][ton1][1]-js1px[x1][i][1])**2)/2
    for ton1 in range(Nma[1]):   
        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
            ypsl11+=Kpn5*((numa[1][ton1][0]+zj2px[x1][i][0]-numa[1][ton1][1]-zj2px[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa[1][ton1][0]+js2px[x1][i][0]-numa[1][ton1][1]-js2px[x1][i][1])**2)/2
    for ton1 in range(Nma1[0]):   
        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
            ypsl11+=Kpn5*((numa1[0][ton1][0]+zj3px[x1][i][0]-numa1[0][ton1][1]-zj3px[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[0][ton1][0]+js3px[x1][i][0]-numa1[0][ton1][1]-js3px[x1][i][1])**2)/2
    for ton1 in range(Nma1[1]):   
        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
            ypsl11+=Kpn5*((numa1[1][ton1][0]+zj4px[x1][i][0]-numa1[1][ton1][1]-zj4px[x1][i][1])**2)/2
            ypsl22+=Kpn5*((numa1[1][ton1][0]+js4px[x1][i][0]-numa1[1][ton1][1]-js4px[x1][i][1])**2)/2 
    return ypsl11,ypsl22


# In[11]:


def energy0():#calculate the energy of the current spindle system
    ypsl0=0
    for x2 in range(pair):
        ypsl0+=Kp3*(((Xkin-l1r[x2])**2+(Xkin1-l4l[x2])**2)/2)+Kp2*((Xpole-l2l[x2])**2+(Xpo-l3r[x2])**2)/2+Kp1*((Xpole-l1l[x2])**2+(Xpo-l4r[x2])**2)/2+Kp4*((Xkin1-Xkin-kdis)**2)/2
        for xn in range(nu[x2]):
            
            if kina[x2][xn]==1 and kinb[x2][xn]==1:
                ppx=motor[x2][xn][0]-motor[x2][xn][1]
                if ppx>xdmotor:
                    ypsl0+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl0+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1[x2]):
            if kinap[x2][xn]+kinbp[x2][xn]==2:
                ppx=motorp[x2][xn][0]-motorp[x2][xn][1]
                if ppx>xdmotor:
                    ypsl0+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl0+=Km*((ppx+xdmotor)**2)/2
        for xn in range(nu1x[x2]):
            if kinapx[x2][xn]+kinbpx[x2][xn]==2:
                ppx=motorpx[x2][xn][0]-motorpx[x2][xn][1]
                if ppx>xdmotor:
                    ypsl0+=Km*((ppx-xdmotor)**2)/2
                if ppx<-xdmotor:
                    ypsl0+=Km*((ppx+xdmotor)**2)/2
    for ton1 in range(Nma[0]):   
        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
            ypsl0+=Kpn5*((numa[0][ton1][0]-numa[0][ton1][1])**2)/2
    for ton1 in range(Nma[1]):   
        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
            ypsl0+=Kpn5*((numa[1][ton1][0]-numa[1][ton1][1])**2)/2
    for ton1 in range(Nma1[0]):   
        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
            ypsl0+=Kpn5*((numa1[0][ton1][0]-numa1[0][ton1][1])**2)/2
    for ton1 in range(Nma1[1]):   
        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
            ypsl0+=Kpn5*((numa1[1][ton1][0]-numa1[1][ton1][1])**2)/2
                
    return ypsl0


# In[12]:


def checkF():#Calculate the force of the two kinetochores, two spindle poles and four MTs

    F0=0#the force of the left spindle pole
    F7=0#the force of the right spindle pole
    F5=0#the force of the left kinetochore
    F6=0#the force of the right kinetochore
    Fcc=0
    for x2 in range(pair):
        F0+=Kp1*(l1l[x2]-Xpole)+Kp2*(l2l[x2]-Xpole)
        F7+=Kp1*(l4r[x2]-Xpo)+Kp2*(l3r[x2]-Xpo)
        F5+=Kp3*(l1r[x2]-Xkin)+Kp4*(Xkin1-Xkin-kdis)
        F6+=Kp4*(Xkin-Xkin1+kdis)+Kp3*(l4l[x2]-Xkin1)
        F3[x2]=Kp2*(Xpo-l3r[x2])#the force of the iMT connected to the right spindle pole

        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)#the force of the kMT connected to the left spindle pole

        F2[x2]=Kp2*(l2l[x2]-Xpole)#the force of the iMT connected to the left spindle pole

        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])#the force of the kMT connected to the right spindle pole

        for ton1 in range(nu[x2]):
            if kina[x2][ton1]==1 and kinb[x2][ton1]==1:
                dxc[x2][ton1]=motor[x2][ton1][0]-motor[x2][ton1][1]
                if dxc[x2][ton1]>xdmotor:
                    F3[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                    F2[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                if dxc[x2][ton1]<-xdmotor:
                    F3[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                    F2[x2]+=(dxc[x2][ton1]+xdmotor)*Km
            if kina[x2][ton1]==0 or kinb[x2][ton1]==0:
                dxc[x2][ton1]=0
        for ton1 in range(nu1[x2]):
            if kinap[x2][ton1]==1 and kinbp[x2][ton1]==1:
                dxcp[x2][ton1]=motorp[x2][ton1][0]-motorp[x2][ton1][1]
                if dxcp[x2][ton1]>xdmotor:
                    F3[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                    F1[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                if dxcp[x2][ton1]<-xdmotor:
                    F3[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                    F1[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
            if kinap[x2][ton1]==0 or kinbp[x2][ton1]==0:
                dxcp[x2][ton1]=0
        for ton1 in range(nu1x[x2]):
            if kinapx[x2][ton1]==1 and kinbpx[x2][ton1]==1:
                dxcpx[x2][ton1]=motorpx[x2][ton1][0]-motorpx[x2][ton1][1]
                if dxcpx[x2][ton1]>xdmotor:
                    F4[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                    F2[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                if dxcpx[x2][ton1]<-xdmotor:
                    F4[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                    F2[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
            if kinapx[x2][ton1]==0 or kinbpx[x2][ton1]==0:
                dxcpx[x2][ton1]=0
    for ton1 in range(Nma[0]):
        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
            F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
            F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
    for ton1 in range(Nma[1]):
        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
            F3[0]-=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
            F3[1]+=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])   
    for ton1 in range(Nma1[0]):
        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
            F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
            F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
    for ton1 in range(Nma1[1]):
        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
            F4[0]-=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
            F4[1]+=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])  
    for x2 in range(pair):
        if abs(F1[x2])+abs(F2[x2])+abs(F3[x2])+abs(F4[x2])>0.0001:
            Fcc=1
    if abs(F0)+abs(F7)+abs(F5)+abs(F6)>0.0001 or Fcc==1:

        de0,de1,de2,de3,de4,de5,de6,de7=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
        return 1,de0,de1,de2,de3,de4,de5,de6,de7
    return 0,0,0,0,0,0,0,0,0


# main


import random
import math
import numpy
import numpy as np
h = 0.001
Km = 0.55# pN/nm
d=8.2
dcc=8.2
Number=2
vp1=6
vpol=vp1/d
xdmotor=d+0.000001
nant=3#[Eg5]
npar=1
pair=2
ka1=0.0004*nant
ka2=0.0004*nant
Kp2=0.1
Kp1=0.1
Kp3=0.1
Kp4=10
lover=6002.4
kdis=1000.4
llap=4001.6
llap1=4001.6
l1l=[-lover for i in range(pair)]#the kMT connected to the left spindle pole
l1r=[-kdis/2 for i in range(pair)]
l2l=[-lover for i in range(pair)]#the iMT connected to the left spindle pole
l2r=[llap ,llap1]
l3l=[-llap ,-llap1]#the iMT connected to the right spindle pole
l3r=[lover for i in range(pair)]
l4l=[kdis/2 for i in range(pair)]#the kMT connected to the right spindle pole
l4r=[lover for i in range(pair)]
Xpole=-lover
Xpo=lover
Xkin=-kdis/2
Xkin1=kdis/2
vp0=vp1/4
Fp0=3.2
#MCAK
kon=0.0001
koff=0.004
konnu=0.0003
Kpn5=0.03
unu=1
kdif=82
vi=5
vk=4
#stai=0.5
#stak=0.5
stad=0.5
##################
Ximt2=[0  for i in range(pair)]
Ximt3=[0  for i in range(pair)]
Xkmt1=[0  for i in range(pair)]
Xkmt4=[0  for i in range(pair)]
desi2=[0  for i in range(pair)]
desi3=[0  for i in range(pair)]
desk1=[0  for i in range(pair)]
desk4=[0  for i in range(pair)]
pos1=[0  for i in range(pair)]
pos4=[0  for i in range(pair)]
pos2=[0  for i in range(pair)]
desia2=[0  for i in range(pair)]
desia3=[0  for i in range(pair)]
deska1=[0  for i in range(pair)]
deska4=[0  for i in range(pair)]
posa1=[0  for i in range(pair)]
posa4=[0  for i in range(pair)]
posa2=[0  for i in range(pair)]
ximta2=[0  for i in range(pair)]
ximta3=[0  for i in range(pair)]
xkmta1=[0  for i in range(pair)]
xkmta4=[0  for i in range(pair)]
ads=1
adeps1=[0 for i in range(pair)]
adeps2=[0 for i in range(pair)]
adeps3=[0 for i in range(pair)]
adeps4=[0 for i in range(pair)]
nu=[0 for i in range(pair)]
nu1=[0 for i in range(pair)]
nu1x=[0 for i in range(pair)]
monum=[0 for i in range(pair)]
monump=[0 for i in range(pair)]
monumpx=[0 for i in range(pair)]
numc1=[0 for i in range(pair)]
numc2=[0 for i in range(pair)]
numc3=[0 for i in range(pair)]
numc4=[0 for i in range(pair)]
lxian=[0 for i in range(pair)]
lxianp=[0 for i in range(pair)]
lxianpx=[0 for i in range(pair)]
F1=[0 for i in range(pair)]
F2=[0 for i in range(pair)]
F3=[0 for i in range(pair)]
F4=[0 for i in range(pair)]
panx=0
pan1=1
site=round(lover/d)-1
parpair=2
kina =[[]for i in range(pair)]
kinb =[[]for i in range(pair)]
kinap =[[]for i in range(pair)]
kinbp =[[]for i in range(pair)]
kinapx=[[]for i in range(pair)]
kinbpx=[[]for i in range(pair)]
kinam=[[]for i in range(parpair)]
kinbm=[[]for i in range(parpair)]
Nma=[0 for i in range(parpair)]
kinam1=[[]for i in range(parpair)]
kinbm1=[[]for i in range(parpair)]
Nma1=[0 for i in range(parpair)]
numa=[[] for i in range(parpair)]
numa1=[[] for i in range(parpair)]
motor =[[]for i in range(pair)]
motorp=[[]for i in range(pair)]
motorpx=[[]for i in range(pair)]
motormc1=[[]for i in range(pair)]
motormc2=[[]for i in range(pair)]
motormc3=[[]for i in range(pair)]
motormc4=[[]for i in range(pair)]
dxc =[[]for i in range(pair)]
dxcp =[[]for i in range(pair)]
dxcpx =[[]for i in range(pair)]
pan=[[]for i in range(pair)]
ypsl1=[[]for i in range(pair)]
ypsl2=[[]for i in range(pair)]
panp=[[]for i in range(pair)]
ypsl1p=[[]for i in range(pair)]
ypsl2p=[[]for i in range(pair)]
panpx=[[]for i in range(pair)]
ypsl1px=[[]for i in range(pair)]
ypsl2px=[[]for i in range(pair)]
zj0=[[] for i in range(pair)]
js0=[[]for i in range(pair)]
zj1=[[]for i in range(pair)]
js1=[[]for i in range(pair)]
zj2=[[]for i in range(pair)]
js2=[[]for i in range(pair)]
zj3=[[]for i in range(pair)]
js3=[[]for i in range(pair)]
zj4=[[]for i in range(pair)]
js4=[[]for i in range(pair)]
zj5=[[]for i in range(pair)]
js5=[[]for i in range(pair)]
zj6=[[]for i in range(pair)]
js6=[[]for i in range(pair)]
zj7=[[]for i in range(pair)]
js7=[[]for i in range(pair)]
zj0p=[[]for i in range(pair)]
js0p=[[]for i in range(pair)]
zj1p=[[]for i in range(pair)]
js1p=[[]for i in range(pair)]
zj2p=[[]for i in range(pair)]
js2p=[[]for i in range(pair)]
zj3p=[[]for i in range(pair)]
js3p=[[]for i in range(pair)]
zj4p=[[]for i in range(pair)]
js4p=[[]for i in range(pair)]
zj5p=[[]for i in range(pair)]
js5p=[[]for i in range(pair)]
zj6p=[[]for i in range(pair)]
js6p=[[]for i in range(pair)]
zj7p=[[]for i in range(pair)]
js7p=[[]for i in range(pair)]
zj0px=[[]for i in range(pair)]
js0px=[[]for i in range(pair)]
zj1px=[[]for i in range(pair)]
js1px=[[]for i in range(pair)]
zj2px=[[]for i in range(pair)]
js2px=[[]for i in range(pair)]
zj3px=[[]for i in range(pair)]
js3px=[[]for i in range(pair)]
zj4px=[[]for i in range(pair)]
js4px=[[]for i in range(pair)]
zj5px=[[]for i in range(pair)]
js5px=[[]for i in range(pair)]
zj6px=[[]for i in range(pair)]
js6px=[[]for i in range(pair)]
zj7px=[[]for i in range(pair)]
js7px=[[]for i in range(pair)]

for x1 in range(pair):
    xx=1/8
    for xn in range(2*7):
        xx+=1/2/8
        cx0=l3r[x1]-round((lover-llap)/d*random.uniform(0,xx))*d####??????
        motormc3[x1].append(cx0)
    xx=1/8
    for xn in range(2*7):
        xx+=1/2/8
        cx0=l1l[x1]+round((lover-llap)/d*random.uniform(0,xx))*d
        motormc1[x1].append(cx0)
    xx=1/8
    for xn in range(2*7):
        xx+=1/2/8
        cx0=l4r[x1]-round((lover-llap)/d*random.uniform(0,xx))*d
        motormc4[x1].append(cx0)
    xx=1/8
    for xn in range(2*7):
        xx+=1/2/8
        cx0=l2l[x1]+round((lover-llap)/d*random.uniform(0,xx))*d
        motormc2[x1].append(cx0)
    for xn in range(round(Number*random.uniform(nant*2,nant*2+1))):
        cx0=l3l[x1]+round(lover/d*random.uniform(0,1))*d
        motor[x1].append([cx0,cx0])
        kina[x1].append(round(random.uniform(0,1)))
        kinb[x1].append(round(random.uniform(0,1)))
        ypsl1[x1].append(0)
        dxc[x1].append(0)
        pan[x1].append(0)
        ypsl2[x1].append(0)
        zj0[x1].append(0)
        js0[x1].append(0)
        zj1[x1].append([0 for i in range(pair)])
        js1[x1].append([0 for i in range(pair)])
        zj2[x1].append([0 for i in range(pair)])
        js2[x1].append([0 for i in range(pair)])
        zj3[x1].append([0 for i in range(pair)])
        js3[x1].append([0 for i in range(pair)])
        zj4[x1].append([0 for i in range(pair)])
        js4[x1].append([0 for i in range(pair)])
        zj5[x1].append(0)
        js5[x1].append(0)
        zj6[x1].append(0)
        js6[x1].append(0)
        zj7[x1].append(0)
        js7[x1].append(0)

    for xn in range(round(Number*random.uniform(nant,nant+1))):
        cx0=l3l[x1]+round(lover/2/d*random.uniform(0,1))*d
        motorp[x1].append([cx0,cx0])
        kinap[x1].append(round(random.uniform(0,1)))
        kinbp[x1].append(round(random.uniform(0,1)))
        dxcp[x1].append(0)
        panp[x1].append(0)
        ypsl1p[x1].append(0)
        ypsl2p[x1].append(0)
        zj0p[x1].append(0)
        js0p[x1].append(0)
        zj1p[x1].append([0 for i in range(pair)])
        js1p[x1].append([0 for i in range(pair)])
        zj2p[x1].append([0 for i in range(pair)])
        js2p[x1].append([0 for i in range(pair)])
        zj3p[x1].append([0 for i in range(pair)])
        js3p[x1].append([0 for i in range(pair)])
        zj4p[x1].append([0 for i in range(pair)])
        js4p[x1].append([0 for i in range(pair)])
        zj5p[x1].append(0)
        js5p[x1].append(0)
        zj6p[x1].append(0)
        js6p[x1].append(0)
        zj7p[x1].append(0)
        js7p[x1].append(0)
    
    for xn in range(round(Number*random.uniform(nant,nant+1))):
        cx0=l4l[x1]+round(lover/2/d*random.uniform(0,1))*d
        motorpx[x1].append([cx0,cx0])
        kinapx[x1].append(round(random.uniform(0,1)))
        kinbpx[x1].append(round(random.uniform(0,1)))
        dxcpx[x1].append(0)
        panpx[x1].append(0)
        ypsl1px[x1].append(0)
        ypsl2px[x1].append(0)
        zj0px[x1].append(0)
        js0px[x1].append(0)
        zj1px[x1].append([0 for i in range(pair)])
        js1px[x1].append([0 for i in range(pair)])
        zj2px[x1].append([0 for i in range(pair)])
        js2px[x1].append([0 for i in range(pair)])
        zj3px[x1].append([0 for i in range(pair)])
        js3px[x1].append([0 for i in range(pair)])
        zj4px[x1].append([0 for i in range(pair)])
        js4px[x1].append([0 for i in range(pair)])
        zj5px[x1].append(0)
        js5px[x1].append(0)
        zj6px[x1].append(0)
        js6px[x1].append(0)
        zj7px[x1].append(0)
        js7px[x1].append(0)

for xn in range(16):
    kinam[0].append(1)
    kinbm[0].append(1)
    cx0=l2l[0]+round(lover/2/d*random.uniform(0,1))*d
    numa[0].append([cx0,cx0])
for xn in range(16):
    kinam[1].append(1)
    kinbm[1].append(1)
    cx0=l3r[0]-round(lover/2/d*random.uniform(0,1))*d
    numa[1].append([cx0,cx0])

Nma[0]=len(kinam[0])
Nma[1]=len(kinam[1])
Nma1[0]=len(kinam1[0])
Nma1[1]=len(kinam1[1])
nunum=[0 for i in range(2)]
nunum1=[0 for i in range(2)]
for xbb in range(Nma[0]):
    if kinam[0][xbb]==1 and kinbm[0][xbb]==1:
        nunum[0]+=1
for xbb in range(Nma[1]):
    if kinam[1][xbb]==1 and kinbm[1][xbb]==1:
        nunum[1]+=1
for xbb in range(Nma1[0]):
    if kinam1[0][xbb]==1 and kinbm1[0][xbb]==1:
        nunum1[0]+=1
for xbb in range(Nma1[1]):
    if kinam1[1][xbb]==1 and kinbm1[1][xbb]==1:
        nunum1[1]+=1
for x1 in range(pair):
    nu[x1]=len(kina[x1])
    nu1[x1]=len(kinap[x1])
    nu1x[x1]=len(kinapx[x1])
ypsl0=energy0()
ypsl0p=ypsl0
ypsl0px=ypsl0
for x1 in range(pair):
    monum[x1]=sum(kina[x1])
    for xn in range(nu[x1]):
        if kina[x1][xn]==1 and kinb[x1][xn]==0:
            monum[x1]-=1
    monump[x1]=sum(kinap[x1])
    for xn in range(nu1[x1]):
        if kinap[x1][xn]==1 and kinbp[x1][xn]==0:
            monump[x1]-=1
    monumpx[x1]=sum(kinapx[x1])
    for xn in range(nu1x[x1]):
        if kinapx[x1][xn]==1 and kinbpx[x1][xn]==0:
            monumpx[x1]-=1


# loop


for a in range(0,20000000):
    
    ran4 = random.uniform(0, 1)
    ss=min(l2r[0]-l2l[0],l2r[1]-l2l[1])
    if ran4 < konnu*round(ss/d-sum(kinam[0]))*h:
        cx0=l2l[0]+round((l2r[0]-l2l[0])*random.uniform(0,1)/d)*d
        numa[0].append([cx0,cx0])
        kinam[0].append(1)
        kinbm[0].append(0)
        Nma[0]+=1
    ran4 = random.uniform(0, 1)
    if ran4 < konnu*round(ss/d-sum(kinbm[0]))*h:
        cx0=l2l[1]+round((l2r[1]-l2l[1])*random.uniform(0,1)/d)*d
        numa[0].append([cx0,cx0])
        kinam[0].append(0)
        kinbm[0].append(1)
        Nma[0]+=1
    ran4 = random.uniform(0, 1)
    ss=min(l3r[0]-l3l[0],l3r[1]-l3l[1])
    if ran4 < konnu*round(ss/d-sum(kinam[1]))*h:
        cx0=l3l[0]+round((l3r[0]-l3l[0])*random.uniform(0,1)/d)*d
        numa[1].append([cx0,cx0])
        kinam[1].append(1)
        kinbm[1].append(0)
        Nma[1]+=1
    ran4 = random.uniform(0, 1)
    if ran4 < konnu*round(ss/d-sum(kinbm[1]))*h:
        cx0=l3l[1]+round((l3r[1]-l3l[1])*random.uniform(0,1)/d)*d
        numa[1].append([cx0,cx0])
        kinam[1].append(0)
        kinbm[1].append(1)
        Nma[1]+=1
    for xbb in range(Nma[0]):
        if xbb>=Nma[0]-0.1:
            continue
        if kinam[0][xbb]==1 and kinbm[0][xbb]==1:
            if numa[0][xbb][0]>l2r[0] or numa[0][xbb][0]<l2l[0] or numa[0][xbb][1]>l2r[1] or numa[0][xbb][1]<l2l[1]:
                kinam[0][xbb]=0;kinbm[0][xbb]=0
                nunum[0]-=1
                panx=1
                pan1=1
                Nma[0]-=1
                del kinam[0][xbb]
                del kinbm[0][xbb]
                del numa[0][xbb]
                continue
        if kinam[0][xbb]==1 and kinbm[0][xbb]==0:
            ran4 = random.uniform(0, 1)
            if ran4 < unu*h:
                lian1=0
                kinbm[0][xbb]=1
                for xn in range(Nma[0]):
                    if kinbm[0][xn]==1 and xn!=xbb and abs(l2l[1]+round((numa[0][xbb][0]-l2l[1])/d)*d-numa[0][xn][1])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l2l[1]+round((numa[0][xbb][0]-l2l[1])/d)*d>=numa[0][xbb][0]:
                        la=(round((numa[0][xbb][0]-l2l[1])/d)-1)*d
                    if l2l[1]+round((numa[0][xbb][0]-l2l[1])/d)*d<numa[0][xbb][0]:
                        la=(round((numa[0][xbb][0]-l2l[1])/d)+1)*d
                    for xn in range(Nma[0]):
                        if kinbm[0][xn]==1 and xn!=xbb and abs(l2l[1]+la-numa[0][xn][1])<1:
                            lian1=2
                            break
                if lian1==0:
                    numa[0][xbb][1]=l2l[1]+round((numa[0][xbb][0]-l2l[1])/d)*d
                if lian1==1:
                    numa[0][xbb][1]=l2l[1]+la
                if lian1==2:
                    kinbm[0][xbb]=0
                if lian1==0 or lian1==1:
                    panx=1
                    nunum[0]+=1
            
        if kinam[0][xbb]==0 and kinbm[0][xbb]==1:
            ran4 = random.uniform(0, 1)
            if ran4 < unu*h:
                lian1=0
                kinam[0][xbb]=1
                for xn in range(Nma[0]):
                    if kinam[0][xn]==1 and xn!=xbb and abs(l2l[0]+round((numa[0][xbb][1]-l2l[0])/d)*d-numa[0][xn][0])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l2l[0]+round((numa[0][xbb][1]-l2l[0])/d)*d>=numa[0][xbb][1]:
                        la=(round((numa[0][xbb][1]-l2l[0])/d)-1)*d
                    if l2l[0]+round((numa[0][xbb][1]-l2l[0])/d)*d<numa[0][xbb][1]:
                        la=(round((numa[0][xbb][1]-l2l[0])/d)+1)*d
                    for xn in range(Nma[0]):
                        if kinbm[0][xn]==1 and xn!=xbb and abs(l2l[0]+la-numa[0][xn][0])<1:
                            lian1=2
                            break
                if lian1==0:
                    numa[0][xbb][0]=l2l[0]+round((numa[0][xbb][1]-l2l[0])/d)*d
                if lian1==1:
                    numa[0][xbb][0]=l2l[0]+la
                if lian1==2:
                    kinam[0][xbb]=0
                if lian1==0 or lian1==1:
                    panx=1
                    nunum[0]+=1
    for xbb in range(Nma[1]):
        if xbb>=Nma[1]-0.1:
            continue
        if kinam[1][xbb]==1 and kinbm[1][xbb]==1:
            if numa[1][xbb][0]>l3r[0] or numa[1][xbb][0]<l3l[0] or numa[1][xbb][1]>l3r[1] or numa[1][xbb][1]<l3l[1]:
                kinam[1][xbb]=0;kinbm[1][xbb]=0
                nunum[1]-=1
                panx=1
                pan1=1
                Nma[1]-=1
                del kinam[1][xbb]
                del kinbm[1][xbb]
                del numa[1][xbb]
                continue
        if kinam[1][xbb]==1 and kinbm[1][xbb]==0:
            ran4 = random.uniform(0, 1)
            if ran4 < unu*h:
                kinbm[1][xbb]=1
                lian1=0
                for xn in range(Nma[1]):
                    if kinbm[1][xn]==1 and xn!=xbb and abs(l3l[1]+round((numa[1][xbb][0]-l3l[1])/d)*d-numa[1][xn][1])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l3l[1]+round((numa[1][xbb][0]-l3l[1])/d)*d>=numa[1][xbb][0]:
                        la=(round((numa[1][xbb][0]-l3l[1])/d)-1)*d
                    if l3l[1]+round((numa[1][xbb][0]-l3l[1])/d)*d<numa[1][xbb][0]:
                        la=(round((numa[1][xbb][0]-l3l[1])/d)+1)*d
                    for xn in range(Nma[1]):
                        if kinbm[1][xn]==1 and xn!=xbb and abs(l3l[1]+la-numa[1][xn][1])<1:
                            lian1=2
                            break
                if lian1==0:
                    numa[1][xbb][1]=l3l[1]+round((numa[1][xbb][0]-l3l[1])/d)*d
                if lian1==1:
                    numa[1][xbb][1]=l3l[1]+la
                if lian1==2:
                    kinbm[1][xbb]=0
                if lian1==0 or lian1==1:
                    panx=1
                    nunum[1]+=1
        if kinam[1][xbb]==0 and kinbm[1][xbb]==1:
            ran4 = random.uniform(0, 1)
            if ran4 < unu*h:
                kinam[1][xbb]=1
                lian1=0
                for xn in range(Nma[1]):
                    if kinam[1][xn]==1 and xn!=xbb and abs(l3l[0]+round((numa[1][xbb][1]-l3l[0])/d)*d-numa[1][xn][0])<1:
                        lian1=1
                        break
                if lian1==1:
                    if l3l[0]+round((numa[1][xbb][1]-l3l[0])/d)*d>=numa[1][xbb][1]:
                        la=(round((numa[1][xbb][1]-l3l[0])/d)-1)*d
                    if l3l[0]+round((numa[1][xbb][1]-l3l[0])/d)*d<numa[1][xbb][1]:
                        la=(round((numa[1][xbb][1]-l3l[0])/d)+1)*d
                    for xn in range(Nma[1]):
                        if kinbm[1][xn]==1 and xn!=xbb and abs(l3l[0]+la-numa[1][xn][0])<1:
                            lian1=2
                            break
                if lian1==0:
                    numa[1][xbb][0]=l3l[0]+round((numa[1][xbb][1]-l3l[0])/d)*d
                if lian1==1:
                    numa[1][xbb][0]=l3l[0]+la
                if lian1==2:
                    kinam[1][xbb]=0
                if lian1==0 or lian1==1:
                    panx=1
                    nunum[1]+=1
        if panx==1:   
        panx=0
        pan1=1
        ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
        if ckF==1:
            for x2 in range(pair):
                for xn in range(nu[x2]):
                    if kina[x2][xn]==1:
                        motor[x2][xn][0]+=de1[x2]
                    if kinb[x2][xn]==1:
                        motor[x2][xn][1]+=de2[x2]
                for xn in range(nu1[x2]):
                    if kinap[x2][xn]==1:
                        motorp[x2][xn][0]+=de3[x2]
                    if kinbp[x2][xn]==1:
                        motorp[x2][xn][1]+=de2[x2]
                for xn in range(nu1x[x2]):
                    if kinapx[x2][xn]==1:
                        motorpx[x2][xn][0]+=de4[x2]
                    if kinbpx[x2][xn]==1:
                        motorpx[x2][xn][1]+=de1[x2]
                for xn in range(numc1[x2]):
                    motormc1[x2][xn]+=de3[x2]
                for xn in range(numc2[x2]):
                    motormc2[x2][xn]+=de1[x2]
                for xn in range(numc3[x2]):
                    motormc3[x2][xn]+=de2[x2]
                for xn in range(numc4[x2]):
                    motormc4[x2][xn]+=de4[x2]
                l2l[x2]+=de1[x2]
                l2r[x2]+=de1[x2]
                l3l[x2]+=de2[x2]
                l3r[x2]+=de2[x2]
                l1l[x2]+=de3[x2]
                l1r[x2]+=de3[x2]
                l4l[x2]+=de4[x2]
                l4r[x2]+=de4[x2]
                Ximt2[x2]+=de1[x2]
                Ximt3[x2]+=de2[x2]
                Xkmt1[x2]+=de3[x2]
                Xkmt4[x2]+=de4[x2]
            for ton1 in range(Nma[0]):   
                if kinam[0][ton1]==1:
                    numa[0][ton1][0]+=de1[0]
                if kinbm[0][ton1]==1:
                    numa[0][ton1][1]+=de1[1]
            for ton1 in range(Nma[1]):   
                if kinam[1][ton1]==1:
                    numa[1][ton1][0]+=de2[0]
                if kinbm[1][ton1]==1:
                    numa[1][ton1][1]+=de2[1]
            for ton1 in range(Nma1[0]):   
                if kinam1[0][ton1]==1:
                    numa1[0][ton1][0]+=de3[0]
                if kinbm1[0][ton1]==1:
                    numa1[0][ton1][1]+=de3[1]
            for ton1 in range(Nma1[1]):   
                if kinam1[1][ton1]==1:
                    numa1[1][ton1][0]+=de4[0]
                if kinbm1[1][ton1]==1:
                    numa1[1][ton1][1]+=de4[1]
            Xpole+=de0
            Xpo+=de7
            Xkin+=de5
            Xkin1+=de6
            ypsl0=energy0()
            ypsl0p=ypsl0
            ypsl0px=ypsl0
            ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
            if ckF==1:
                for x2 in range(pair):
                    for xn in range(nu[x2]):
                        if kina[x2][xn]==1:
                            motor[x2][xn][0]+=de1[x2]
                        if kinb[x2][xn]==1:
                            motor[x2][xn][1]+=de2[x2]
                    for xn in range(nu1[x2]):
                        if kinap[x2][xn]==1:
                            motorp[x2][xn][0]+=de3[x2]
                        if kinbp[x2][xn]==1:
                            motorp[x2][xn][1]+=de2[x2]
                    for xn in range(nu1x[x2]):
                        if kinapx[x2][xn]==1:
                            motorpx[x2][xn][0]+=de4[x2]
                        if kinbpx[x2][xn]==1:
                            motorpx[x2][xn][1]+=de1[x2]
                    for xn in range(numc1[x2]):
                        motormc1[x2][xn]+=de3[x2]
                    for xn in range(numc2[x2]):
                        motormc2[x2][xn]+=de1[x2]
                    for xn in range(numc3[x2]):
                        motormc3[x2][xn]+=de2[x2]
                    for xn in range(numc4[x2]):
                        motormc4[x2][xn]+=de4[x2]
                    l2l[x2]+=de1[x2]
                    l2r[x2]+=de1[x2]
                    l3l[x2]+=de2[x2]
                    l3r[x2]+=de2[x2]
                    l1l[x2]+=de3[x2]
                    l1r[x2]+=de3[x2]
                    l4l[x2]+=de4[x2]
                    l4r[x2]+=de4[x2]
                    Ximt2[x2]+=de1[x2]
                    Ximt3[x2]+=de2[x2]
                    Xkmt1[x2]+=de3[x2]
                    Xkmt4[x2]+=de4[x2]
                for ton1 in range(Nma[0]):   
                    if kinam[0][ton1]==1:
                        numa[0][ton1][0]+=de1[0]
                    if kinbm[0][ton1]==1:
                        numa[0][ton1][1]+=de1[1]
                for ton1 in range(Nma[1]):   
                    if kinam[1][ton1]==1:
                        numa[1][ton1][0]+=de2[0]
                    if kinbm[1][ton1]==1:
                        numa[1][ton1][1]+=de2[1]
                for ton1 in range(Nma1[0]):   
                    if kinam1[0][ton1]==1:
                        numa1[0][ton1][0]+=de3[0]
                    if kinbm1[0][ton1]==1:
                        numa1[0][ton1][1]+=de3[1]
                for ton1 in range(Nma1[1]):   
                    if kinam1[1][ton1]==1:
                        numa1[1][ton1][0]+=de4[0]
                    if kinbm1[1][ton1]==1:
                        numa1[1][ton1][1]+=de4[1]
                Xpole+=de0
                Xpo+=de7
                Xkin+=de5
                Xkin1+=de6
                ypsl0=energy0()
                ypsl0p=ypsl0
                ypsl0px=ypsl0
##########################################################
    for x1 in range(pair):
        lspa1=min(l3l[x1]-l1l[x1],Xkin-l1l[x1])##1MT
        lspa2=l3l[x1]-l2l[x1]##2MT
        lspa3=l3r[x1]-l2r[x1]#3MT
        lspa4=min(l4r[x1]-l2r[x1],l4r[x1]-Xkin1)
        numc1[x1]=len(motormc1[x1])
        numc2[x1]=len(motormc2[x1])
        numc3[x1]=len(motormc3[x1])
        numc4[x1]=len(motormc4[x1])
        nn1=0
        nn2=0
        nn3=0
        for xn in range(numc1[x1]):
            if motormc1[x1][xn]<l1l[x1]+1000:
                nn1+=1
            if motormc1[x1][xn]>=l1l[x1]+1000 and motormc1[x1][xn]<l1l[x1]+2000:
                nn2+=1
            if motormc1[x1][xn]>=l1l[x1]+2000 and motormc1[x1][xn]<l1l[x1]+3000:
                nn3+=1
        if lspa1>1000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn1*d)*h:
                cx0=l1l[x1]+round(1000/d*random.uniform(0,1))*d
                it=0
                while it <numc1[x1]:
                    if abs(cx0-motormc1[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc1[x1].append(cx0) 
        if lspa1>2000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn2*d)*h:
                cx0=l1l[x1]+round(2000/d*random.uniform(0.5,1))*d
                it=0
                while it <numc1[x1]:
                    if abs(cx0-motormc1[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc1[x1].append(cx0) 
        if lspa1>3000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn3*d)*h:
                cx0=l1l[x1]+round(3000/d*random.uniform(0.67,1))*d
                it=0
                while it <numc1[x1]:
                    if abs(cx0-motormc1[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc1[x1].append(cx0)
        if lspa1<=1000:
            nc1=numc1[x1]
            ld=lspa1
        if lspa1<=2000 and lspa1>1000:
            nc1=numc1[x1]-nn1
            ld=lspa1-1000
        if lspa1<=3000 and lspa1>2000:
            nc1=numc1[x1]-nn1-nn2
            ld=lspa1-2000
        if lspa1>3000:
            nc1=numc1[x1]-nn1-nn2-nn3
            ld=lspa1-3000
        ran5 = random.uniform(0, 1)
        if ran5<kon*(ld-nc1*d)*h:
            cx0=l1l[x1]+round(lspa1/d*random.uniform(1-ld/lspa1,1))*d
            it=0
            while it <numc1[x1]:
                if abs(cx0-motormc1[x1][it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc1[x1].append(cx0)
        nn1=0
        nn2=0
        nn3=0
        for xn in range(numc2[x1]):
            if motormc2[x1][xn]<l2l[x1]+1000:
                nn1+=1
            if motormc2[x1][xn]>=l2l[x1]+1000 and motormc2[x1][xn]<l2l[x1]+2000:
                nn2+=1
            if motormc2[x1][xn]>=l2l[x1]+2000 and motormc2[x1][xn]<l2l[x1]+3000:
                nn3+=1
        if lspa2>1000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn1*d)*h:
                cx0=l2l[x1]+round(1000/d*random.uniform(0,1))*d
                it=0
                while it <numc2[x1]:
                    if abs(cx0-motormc2[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc2[x1].append(cx0) 
        if lspa2>2000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn2*d)*h:
                cx0=l2l[x1]+round(2000/d*random.uniform(0.5,1))*d
                it=0
                while it <numc2[x1]:
                    if abs(cx0-motormc2[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc2[x1].append(cx0) 
        if lspa2>3000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn3*d)*h:
                cx0=l2l[x1]+round(3000/d*random.uniform(0.67,1))*d
                it=0
                while it <numc2[x1]:
                    if abs(cx0-motormc2[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc2[x1].append(cx0)
        if lspa2<=1000:
            nc1=numc2[x1]
            ld=lspa2
        if lspa2<=2000 and lspa2>1000:
            nc1=numc2[x1]-nn1
            ld=lspa2-1000
        if lspa2<=3000 and lspa2>2000:
            nc1=numc2[x1]-nn1-nn2
            ld=lspa2-2000
        if lspa2>3000:
            nc1=numc2[x1]-nn1-nn2-nn3
            ld=lspa2-3000
        ran5 = random.uniform(0, 1)
        if ran5<kon*(ld-nc1*d)*h:
            cx0=l2l[x1]+round(lspa2/d*random.uniform(1-ld/lspa2,1))*d
            it=0
            while it <numc2[x1]:
                if abs(cx0-motormc2[x1][it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc2[x1].append(cx0)
        nn1=0
        nn2=0
        nn3=0
        for xn in range(numc3[x1]):
            if motormc3[x1][xn]>l3r[x1]-1000:
                nn1+=1
            if motormc3[x1][xn]<=l3r[x1]-1000 and motormc3[x1][xn]>l3r[x1]-2000:
                nn2+=1
            if motormc3[x1][xn]<=l3r[x1]-2000 and motormc3[x1][xn]>l3r[x1]-3000:
                nn3+=1
        if lspa3>1000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn1*d)*h:
                cx0=l3r[x1]-round(1000/d*random.uniform(0,1))*d
                it=0
                while it <numc3[x1]:
                    if abs(cx0-motormc3[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc3[x1].append(cx0) 
        if lspa3>2000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn2*d)*h:
                cx0=l3r[x1]-round(2000/d*random.uniform(0.5,1))*d
                it=0
                while it <numc3[x1]:
                    if abs(cx0-motormc3[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc3[x1].append(cx0) 
        if lspa3>3000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn3*d)*h:
                cx0=l3r[x1]-round(3000/d*random.uniform(0.67,1))*d
                it=0
                while it <numc3[x1]:
                    if abs(cx0-motormc3[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc3[x1].append(cx0)
        if lspa3<=1000:
            nc1=numc3[x1]
            ld=lspa3
        if lspa3<=2000 and lspa3>1000:
            nc1=numc3[x1]-nn1
            ld=lspa3-1000
        if lspa3<=3000 and lspa3>2000:
            nc1=numc3[x1]-nn1-nn2
            ld=lspa3-2000
        if lspa3>3000:
            nc1=numc3[x1]-nn1-nn2-nn3
            ld=lspa3-3000
        ran4 = random.uniform(0, 1)
        if ran4<kon*(ld-nc1*d)*h:
            cx0=l3r[x1]-round(lspa3/d*random.uniform(1-ld/lspa3,1))*d
            it=0
            while it <numc3[x1]:
                if abs(cx0-motormc3[x1][it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc3[x1].append(cx0)
        nn1=0
        nn2=0
        nn3=0
        for xn in range(numc4[x1]):
            if motormc4[x1][xn]>l4r[x1]-1000:
                nn1+=1
            if motormc4[x1][xn]<=l4r[x1]-1000 and motormc4[x1][xn]>l4r[x1]-2000:
                nn2+=1
            if motormc4[x1][xn]<=l4r[x1]-2000 and motormc4[x1][xn]>l4r[x1]-3000:
                nn3+=1
        if lspa4>1000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn1*d)*h:
                cx0=l4r[x1]-round(1000/d*random.uniform(0,1))*d
                it=0
                while it <numc4[x1]:
                    if abs(cx0-motormc4[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc4[x1].append(cx0) 
        if lspa4>2000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn2*d)*h:
                cx0=l4r[x1]-round(2000/d*random.uniform(0.5,1))*d
                it=0
                while it <numc4[x1]:
                    if abs(cx0-motormc4[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc4[x1].append(cx0) 
        if lspa4>3000:
            ran5 = random.uniform(0, 1)
            if ran5<kon*(1000-nn3*d)*h:
                cx0=l4r[x1]-round(3000/d*random.uniform(0.67,1))*d
                it=0
                while it <numc4[x1]:
                    if abs(cx0-motormc4[x1][it])<0:
                        cx0+=d
                        it=-1
                    it+=1
                motormc4[x1].append(cx0)
        if lspa4<=1000:
            nc1=numc4[x1]
            ld=lspa4
        if lspa4<=2000 and lspa4>1000:
            nc1=numc4[x1]-nn1
            ld=lspa4-1000
        if lspa4<=3000 and lspa4>2000:
            nc1=numc4[x1]-nn1-nn2
            ld=lspa4-2000
        if lspa4>3000:
            nc1=numc4[x1]-nn1-nn2-nn3
            ld=lspa4-3000
        ran4 = random.uniform(0, 1)
        if ran4<kon*(ld-nc1*d)*h:
            cx0=l4r[x1]-round(lspa4/d*random.uniform(1-ld/lspa4,1))*d
            it=0
            while it <numc4[x1]:
                if abs(cx0-motormc4[x1][it])<0:
                    cx0+=d
                    it=-1
                it+=1
            motormc4[x1].append(cx0)
        
##############################################################
        lxian[x1]=l2r[x1]-l3l[x1]
        lxianp[x1]=l1r[x1]-l3l[x1]
        lxianpx[x1]=l2r[x1]-l4l[x1]
        lspa=l2r[x1]-l3l[x1]
        ran4 = random.uniform(0, 1)
        ran5 = random.uniform(0, 1)
        if ran4 < ka1* (round(lspa/d)-1-sum(kina[x1]))* h:
            cx0=l2r[x1]-(round(lspa*random.uniform(0,1)/d)-1)*d
            motor[x1].append([cx0,cx0])
            kina[x1].append(1)
            kinb[x1].append(0)
            ypsl1[x1].append(0)
            dxc[x1].append(0)
            ypsl2[x1].append(0)
            zj0[x1].append(0)
            js0[x1].append(0)
            zj1[x1].append([0 for i in range(pair)])
            js1[x1].append([0 for i in range(pair)])
            zj2[x1].append([0 for i in range(pair)])
            js2[x1].append([0 for i in range(pair)])
            zj3[x1].append([0 for i in range(pair)])
            js3[x1].append([0 for i in range(pair)])
            zj4[x1].append([0 for i in range(pair)])
            js4[x1].append([0 for i in range(pair)])
            zj5[x1].append(0)
            js5[x1].append(0)
            zj6[x1].append(0)
            js6[x1].append(0)
            zj7[x1].append(0)
            js7[x1].append(0)
            pan[x1].append(0)
            nu[x1]+=1
        if ran5 < ka1* (round(lspa/d)-1-sum(kinb[x1]))* h:
            dxc[x1].append(0)
            kina[x1].append(0)
            kinb[x1].append(1)
            cx0=l3l[x1]+(round(lspa*random.uniform(0,1)/d)-1)*d
            motor[x1].append([cx0,cx0])
            ypsl1[x1].append(0)
            ypsl2[x1].append(0)
            zj0[x1].append(0)
            js0[x1].append(0)
            zj1[x1].append([0 for i in range(pair)])
            js1[x1].append([0 for i in range(pair)])
            zj2[x1].append([0 for i in range(pair)])
            js2[x1].append([0 for i in range(pair)])
            zj3[x1].append([0 for i in range(pair)])
            js3[x1].append([0 for i in range(pair)])
            zj4[x1].append([0 for i in range(pair)])
            js4[x1].append([0 for i in range(pair)])
            zj5[x1].append(0)
            js5[x1].append(0)
            zj6[x1].append(0)
            js6[x1].append(0)
            zj7[x1].append(0)
            js7[x1].append(0)
            pan[x1].append(0)
            nu[x1]+=1
        ran4 = random.uniform(0, 1)
        ran5 = random.uniform(0, 1)
        if ran4 < ka2* (round(lxianp[x1]/d)-1-sum(kinap[x1]))* h:
            cx0=l1r[x1]-(round(lxianp[x1]*random.uniform(0,1)/d)-1)*d
            motorp[x1].append([cx0,cx0])
            kinap[x1].append(1)
            kinbp[x1].append(0)
            dxcp[x1].append(0)
            ypsl1p[x1].append(0)
            ypsl2p[x1].append(0)
            zj0p[x1].append(0)
            js0p[x1].append(0)
            zj1p[x1].append([0 for i in range(pair)])
            js1p[x1].append([0 for i in range(pair)])
            zj2p[x1].append([0 for i in range(pair)])
            js2p[x1].append([0 for i in range(pair)])
            zj3p[x1].append([0 for i in range(pair)])
            js3p[x1].append([0 for i in range(pair)])
            zj4p[x1].append([0 for i in range(pair)])
            js4p[x1].append([0 for i in range(pair)])
            zj5p[x1].append(0)
            js5p[x1].append(0)
            zj6p[x1].append(0)
            js6p[x1].append(0)
            zj7p[x1].append(0)
            js7p[x1].append(0)
            panp[x1].append(0)
            nu1[x1]+=1
        if ran5 < ka2* (round(lxianp[x1]/d)-1-sum(kinbp[x1]))* h:
            cx0=l3l[x1]+(round(lxianp[x1]*random.uniform(0,1)/d)-1)*d
            motorp[x1].append([cx0,cx0])
            kinap[x1].append(0)
            kinbp[x1].append(1)
            dxcp[x1].append(0)
            ypsl1p[x1].append(0)
            ypsl2p[x1].append(0)               
            panp[x1].append(0)
            zj0p[x1].append(0)
            js0p[x1].append(0)
            zj1p[x1].append([0 for i in range(pair)])
            js1p[x1].append([0 for i in range(pair)])
            zj2p[x1].append([0 for i in range(pair)])
            js2p[x1].append([0 for i in range(pair)])
            zj3p[x1].append([0 for i in range(pair)])
            js3p[x1].append([0 for i in range(pair)])
            zj4p[x1].append([0 for i in range(pair)])
            js4p[x1].append([0 for i in range(pair)])
            zj5p[x1].append(0)
            js5p[x1].append(0)
            zj6p[x1].append(0)
            js6p[x1].append(0)
            zj7p[x1].append(0)
            js7p[x1].append(0)
            nu1[x1]+=1
        ran4 = random.uniform(0, 1)
        ran5 = random.uniform(0, 1)
        if ran4 < ka2* (round(lxianpx[x1]/d)-1-sum(kinapx[x1]))* h:
            cx0=l4l[x1]+(round(lxianpx[x1]*random.uniform(0,1)/d)-1)*d
            motorpx[x1].append([cx0,cx0])
            kinapx[x1].append(1)
            kinbpx[x1].append(0)
            dxcpx[x1].append(0)
            ypsl1px[x1].append(0)
            ypsl2px[x1].append(0)
            zj0px[x1].append(0)
            js0px[x1].append(0)
            zj1px[x1].append([0 for i in range(pair)])
            js1px[x1].append([0 for i in range(pair)])
            zj2px[x1].append([0 for i in range(pair)])
            js2px[x1].append([0 for i in range(pair)])
            zj3px[x1].append([0 for i in range(pair)])
            js3px[x1].append([0 for i in range(pair)])
            zj4px[x1].append([0 for i in range(pair)])
            js4px[x1].append([0 for i in range(pair)])
            zj5px[x1].append(0)
            js5px[x1].append(0)
            zj6px[x1].append(0)
            js6px[x1].append(0)
            zj7px[x1].append(0)
            js7px[x1].append(0)
            panpx[x1].append(0)
            nu1x[x1]+=1
        if ran5 < ka2* (round(lxianpx[x1]/d)-1-sum(kinbpx[x1]))* h:
            cx0=l2r[x1]-(round(lxianpx[x1]*random.uniform(0,1)/d)-1)*d
            motorpx[x1].append([cx0,cx0])
            kinapx[x1].append(0)
            kinbpx[x1].append(1)
            dxcpx[x1].append(0)
            ypsl1px[x1].append(0)
            ypsl2px[x1].append(0)
            panpx[x1].append(0)
            zj0px[x1].append(0)
            js0px[x1].append(0)
            zj1px[x1].append([0 for i in range(pair)])
            js1px[x1].append([0 for i in range(pair)])
            zj2px[x1].append([0 for i in range(pair)])
            js2px[x1].append([0 for i in range(pair)])
            zj3px[x1].append([0 for i in range(pair)])
            js3px[x1].append([0 for i in range(pair)])
            zj4px[x1].append([0 for i in range(pair)])
            js4px[x1].append([0 for i in range(pair)])
            zj5px[x1].append(0)
            js5px[x1].append(0)
            zj6px[x1].append(0)
            js6px[x1].append(0)
            zj7px[x1].append(0)
            js7px[x1].append(0)
            nu1x[x1]+=1
###############################################################    
    for x1 in range(pair):
        ran4 = random.uniform(0, 1)
        ran5 = random.uniform(0, 1)
        numc1[x1]=len(motormc1[x1])
        numc2[x1]=len(motormc2[x1])
        numc3[x1]=len(motormc3[x1])
        numc4[x1]=len(motormc4[x1])
        for i in range(numc3[x1]):
            ran4=random.uniform(0, 1)
            if adeps3[x1]!=0 and motormc3[x1].index(max(motormc3[x1]))==i:
                continue
            if ran4<kdif*h:
                zq=0
                if motormc3[x1][i]+d>=motormc3[x1][motormc3[x1].index(max(motormc3[x1]))]-0.1 and  motormc3[x1][i]<motormc3[x1][motormc3[x1].index(max(motormc3[x1]))] and adeps3[x1]!=0:
                    zq=3
                for xn in range(numc3[x1]):
                    if motormc3[x1][i]+d>=motormc3[x1][xn]-0.1 and  motormc3[x1][i]<motormc3[x1][xn]:
                        zq=3
                if zq==0:
                    motormc3[x1][i]+=d
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc3[x1][i]-d<l2r[x1]:
                    zq=3
                for xn in range(numc3[x1]):
                    if motormc3[x1][i]-d<=motormc3[x1][xn]+0.1 and  motormc3[x1][i]>motormc3[x1][xn]:
                        zq=3
                if zq==0:
                    motormc3[x1][i]-=d
        for i in range(numc1[x1]):
            if adeps1[x1]!=0 and motormc1[x1].index(min(motormc1[x1]))==i:
                continue
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc1[x1][i]+d>l3l[x1] or motormc1[x1][i]+d>Xkin:
                    zq=3
                for xn in range(numc1[x1]):
                    if motormc1[x1][i]+d>=motormc1[x1][xn]-0.1 and  motormc1[x1][i]<motormc1[x1][xn]:
                        zq=3
                if zq==0:
                    motormc1[x1][i]+=d
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc1[x1][i]-d<=motormc1[x1][motormc1[x1].index(min(motormc1[x1]))]+0.1 and  motormc1[x1][i]>motormc1[x1][motormc1[x1].index(min(motormc1[x1]))] and adeps1[x1]!=0:
                    zq=3
                for xn in range(numc1[x1]):
                    if motormc1[x1][i]-d<=motormc1[x1][xn]+0.1 and  motormc1[x1][i]>motormc1[x1][xn]:
                        zq=3
                if zq==0:
                    motormc1[x1][i]-=d
        for i in range(numc4[x1]):
            ran4=random.uniform(0, 1)
            if adeps4[x1]!=0 and motormc4[x1].index(max(motormc4[x1]))==i:
                continue
            if ran4<kdif*h:
                zq=0
                if motormc4[x1][i]+d>=motormc4[x1][motormc4[x1].index(max(motormc4[x1]))]-0.1 and  motormc4[x1][i]<motormc4[x1][motormc4[x1].index(max(motormc4[x1]))] and adeps4[x1]!=0:
                    zq=3
                for xn in range(numc4[x1]):
                    if motormc4[x1][i]+d>=motormc4[x1][xn]-0.1 and  motormc4[x1][i]<motormc4[x1][xn]:
                        zq=3
                if zq==0:
                    motormc4[x1][i]+=d
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc4[x1][i]-d<l2r[x1] or motormc4[x1][i]-d<Xkin1:
                    zq=3
                for xn in range(numc4[x1]):
                    if motormc4[x1][i]-d<=motormc4[x1][xn]+0.1 and  motormc4[x1][i]>motormc4[x1][xn]:
                        zq=3
                if zq==0:
                    motormc4[x1][i]-=d
        for i in range(numc2[x1]):
            if adeps2[x1]!=0 and motormc2[x1].index(min(motormc2[x1]))==i:
                continue
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc2[x1][i]+d>l3l[x1]:
                    zq=3
                for xn in range(numc2[x1]):
                    if motormc2[x1][i]+d>=motormc2[x1][xn]-0.1 and  motormc2[x1][i]<motormc2[x1][xn]:
                        zq=3
                if zq==0:
                    motormc2[x1][i]+=d
            ran4=random.uniform(0, 1)
            if ran4<kdif*h:
                zq=0
                if motormc2[x1][i]-d<=motormc2[x1][motormc2[x1].index(min(motormc2[x1]))]+0.1 and  motormc2[x1][i]>motormc2[x1][motormc2[x1].index(min(motormc2[x1]))] and adeps2[x1]!=0:
                    zq=3
                for xn in range(numc2[x1]):
                    if motormc2[x1][i]-d<=motormc2[x1][xn]+0.1 and  motormc2[x1][i]>motormc2[x1][xn]:
                        zq=3
                if zq==0:
                    motormc2[x1][i]-=d
        if numc3[x1]>0: 
            if adeps3[x1]!=0 and max(motormc3[x1])+3/2*d<l3r[x1]:
                adeps3[x1]=0
            if max(motormc3[x1])+3/2*d>=l3r[x1]:
                adeps3[x1]+=1
                ran4=random.uniform(0, 1)
                if l3r[x1]>=Xpo:
                    vd=vi
                    #stad=stai
                if l3r[x1]<Xpo:
                    vd=vk
                    #stad=stak
                if ran4<vd*h:
                    motormc3[x1][motormc3[x1].index(max(motormc3[x1]))]-=d
                    l3r[x1]-=d
                    desi3[x1]+=d
                    pan1=1
                    panx=1
                ran4=random.uniform(0, 1)
                if ran4<stad*h:
                    motormc3[x1].pop(motormc3[x1].index(max(motormc3[x1])))
                    adeps3[x1]=0
                    numc3[x1]-=1
        if numc1[x1]>0:
            if adeps1[x1]!=0 and min(motormc1[x1])-3/2*d>l1l[x1]:
                adeps1[x1]=0
            if min(motormc1[x1])-3/2*d<=l1l[x1]:
                adeps1[x1]+=1
                ran4=random.uniform(0, 1)
                if l1l[x1]<=Xpole:
                    vd=vi
                    #stad=stai
                if l1l[x1]>Xpole:
                    vd=vk
                    #stad=stak
                if ran4<vd*h:
                    motormc1[x1][motormc1[x1].index(min(motormc1[x1]))]+=d
                    pan1=1
                    l1l[x1]+=d
                    desk1[x1]+=d
                    panx=1
                ran4=random.uniform(0, 1)
                if ran4<stad*h:
                    motormc1[x1].pop(motormc1[x1].index(min(motormc1[x1])))
                    adeps1[x1]=0
                    numc1[x1]-=1
        if numc4[x1]>0: 
            if adeps4[x1]!=0 and max(motormc4[x1])+3/2*d<l4r[x1]:
                adeps4[x1]=0
            if max(motormc4[x1])+3/2*d>=l4r[x1]:
                adeps4[x1]+=1
                ran4=random.uniform(0, 1)
                if l4r[x1]>=Xpo:
                    vd=vi
                    #stad=stai
                if l4r[x1]<Xpo:
                    vd=vk
                    #stad=stak
                if ran4<vd*h:
                    motormc4[x1][motormc4[x1].index(max(motormc4[x1]))]-=d
                    l4r[x1]-=d
                    desk4[x1]+=d
                    pan1=1
                    panx=1
                ran4=random.uniform(0, 1)
                if ran4<stad*h:
                    motormc4[x1].pop(motormc4[x1].index(max(motormc4[x1])))
                    adeps4[x1]=0
                    numc4[x1]-=1
        if numc2[x1]>0:
            if adeps2[x1]!=0 and min(motormc2[x1])-3/2*d>l2l[x1]:
                adeps2[x1]=0
            if min(motormc2[x1])-3/2*d<=l2l[x1]:
                adeps2[x1]+=1
                ran4=random.uniform(0, 1)
                if l2l[x1]<=Xpole:
                    vd=vi
                    #stad=stai
                if l2l[x1]>Xpole:
                    vd=vk
                    #stad=stak
                if ran4<vd*h:
                    motormc2[x1][motormc2[x1].index(min(motormc2[x1]))]+=d
                    pan1=1
                    l2l[x1]+=d
                    desi2[x1]+=d
                    panx=1
                ran4=random.uniform(0, 1)
                if ran4<stad*h:
                    motormc2[x1].pop(motormc2[x1].index(min(motormc2[x1])))
                    adeps2[x1]=0
                    numc2[x1]-=1 

        ran5 = random.uniform(0, 1)
        vpolkmt=vp0*(1+(Kp3*(Xkin-l1r[x1])/Fp0))/d
        if ran5<abs(vpolkmt)*h:
            if vpolkmt>0:
                l1r[x1]+=d
                pos1[x1]+=d
            if vpolkmt<0:
                l1r[x1]-=d
                pos1[x1]-=d
            pan1=1
            panx=1
        ran5 = random.uniform(0, 1)
        #vpolkmt=vp0*(1+(Kp3*(l4l[x1]-Xkin1)/Fp0))/d
        if ran5<abs(vpolkmt)*h:
            if vpolkmt>0:
                l4l[x1]-=d
                pos4[x1]-=d
            if vpolkmt<0:
                l4l[x1]+=d
                pos4[x1]+=d
            pan1=1
            panx=1
        nu[x1]=len(kina[x1])
        nu1[x1]=len(kinap[x1])
        nu1x[x1]=len(kinapx[x1])
        if pan1==1 and panx==1:
            #print('mtd')
            panx=0
            ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
            for x2 in range(pair):
                for xn in range(nu[x2]):
                    if kina[x2][xn]==1:
                        motor[x2][xn][0]+=de1[x2]
                    if kinb[x2][xn]==1:
                        motor[x2][xn][1]+=de2[x2]
                for xn in range(nu1[x2]):
                    if kinap[x2][xn]==1:
                        motorp[x2][xn][0]+=de3[x2]
                    if kinbp[x2][xn]==1:
                        motorp[x2][xn][1]+=de2[x2]
                for xn in range(nu1x[x2]):
                    if kinapx[x2][xn]==1:
                        motorpx[x2][xn][0]+=de4[x2]
                    if kinbpx[x2][xn]==1:
                        motorpx[x2][xn][1]+=de1[x2]
                for xn in range(numc1[x2]):
                    motormc1[x2][xn]+=de3[x2]
                for xn in range(numc2[x2]):
                    motormc2[x2][xn]+=de1[x2]
                for xn in range(numc3[x2]):
                    motormc3[x2][xn]+=de2[x2]
                for xn in range(numc4[x2]):
                    motormc4[x2][xn]+=de4[x2]
                l2l[x2]+=de1[x2]
                l2r[x2]+=de1[x2]
                l3l[x2]+=de2[x2]
                l3r[x2]+=de2[x2]
                l1l[x2]+=de3[x2]
                l1r[x2]+=de3[x2]
                l4l[x2]+=de4[x2]
                l4r[x2]+=de4[x2]
                Ximt2[x2]+=de1[x2]
                Ximt3[x2]+=de2[x2]
                Xkmt1[x2]+=de3[x2]
                Xkmt4[x2]+=de4[x2]
            for ton1 in range(Nma[0]):   
                if kinam[0][ton1]==1:
                    numa[0][ton1][0]+=de1[0]
                if kinbm[0][ton1]==1:
                    numa[0][ton1][1]+=de1[1]
            for ton1 in range(Nma[1]):   
                if kinam[1][ton1]==1:
                    numa[1][ton1][0]+=de2[0]
                if kinbm[1][ton1]==1:
                    numa[1][ton1][1]+=de2[1]
            for ton1 in range(Nma1[0]):   
                if kinam1[0][ton1]==1:
                    numa1[0][ton1][0]+=de3[0]
                if kinbm1[0][ton1]==1:
                    numa1[0][ton1][1]+=de3[1]
            for ton1 in range(Nma1[1]):   
                if kinam1[1][ton1]==1:
                    numa1[1][ton1][0]+=de4[0]
                if kinbm1[1][ton1]==1:
                    numa1[1][ton1][1]+=de4[1]
            Xpole+=de0
            Xpo+=de7
            Xkin+=de5
            Xkin1+=de6
            ypsl0=energy0()
            ypsl0p=ypsl0
            ypsl0px=ypsl0
            
    #################MT2, MT3, polymerization
    for x1 in range(pair):
        ran4=random.uniform(0, 1)
        if ran4<vpol*h:
            l2r[x1]+=d
            pos2[x1]+=d
            if l2r[x1]>l3r[x1] or l2r[x1]>Xpo:
                l2r[x1]-=d
                pos2[x1]-=d
        ran4=random.uniform(0, 1)
        if ran4<vpol*h:
            l3l[x1]-=d
            if l3l[x1]<l2l[x1] or l3l[x1]<Xpole:
                l3l[x1]+=d
# MCAK dissociation
        i=0
        while i < numc1[x1]:
            if adeps1[x1]!=0 and motormc1[x1].index(min(motormc1[x1]))==i:
                i=i+1
                continue
            ran4=random.uniform(0, 1)
            if ran4<koff*h:
                del motormc1[x1][i]
                numc1[x1]-=1
                continue
            if motormc1[x1][i]<l1l[x1] or motormc1[x1][i]>l1r[x1]:
                del motormc1[x1][i]
                numc1[x1]-=1
                i-=1
            i+=1
        i=0
        while i < numc3[x1]:
            if adeps3[x1]!=0 and motormc3[x1].index(max(motormc3[x1]))==i:
                i=i+1
                continue
            ran4=random.uniform(0, 1)
            if ran4<koff*h:
                del motormc3[x1][i]
                numc3[x1]-=1
                continue
            if motormc3[x1][i]<l3l[x1] or motormc3[x1][i]>l3r[x1]:
                del motormc3[x1][i]
                numc3[x1]-=1
                i-=1
            i+=1

        i=0
        while i < numc2[x1]:
            if adeps2[x1]!=0 and motormc2[x1].index(min(motormc2[x1]))==i:
                i=i+1
                continue
            ran4=random.uniform(0, 1)
            if ran4<koff*h:
                del motormc2[x1][i]
                numc2[x1]-=1
                continue
            if motormc2[x1][i]<l2l[x1] or motormc2[x1][i]>l2r[x1]:
                del motormc2[x1][i]
                numc2[x1]-=1
                i-=1
            i+=1
        i=0
        while i < numc4[x1]:
            if adeps4[x1]!=0 and motormc4[x1].index(max(motormc4[x1]))==i:
                i=i+1
                continue
            ran4=random.uniform(0, 1)
            if ran4<koff*h:
                del motormc4[x1][i]
                numc4[x1]-=1
                continue
            if motormc4[x1][i]<l4l[x1] or motormc4[x1][i]>l4r[x1]:
                del motormc4[x1][i]
                numc4[x1]-=1
                i-=1    
            i+=1
####################################################
    if pan1==1:
        for x1 in range(pair):
            monum[x1]=sum(kina[x1])
            for xn in range(nu[x1]):
                if kina[x1][xn]==1 and kinb[x1][xn]==0:
                    monum[x1]-=1
            monump[x1]=sum(kinap[x1])
            for xn in range(nu1[x1]):
                if kinap[x1][xn]==1 and kinbp[x1][xn]==0:
                    monump[x1]-=1
            monumpx[x1]=sum(kinapx[x1])
            for xn in range(nu1x[x1]):
                if kinapx[x1][xn]==1 and kinbpx[x1][xn]==0:
                    monumpx[x1]-=1
    for x1 in range(pair):
        for i in range(0,nu[x1]):
            if kina[x1][i]+kinb[x1][i]<2:
                if pan1==1:
                    pan[x1][i]=0
                if kina[x1][i]==1 and kinb[x1][i]==0:
                    na=stepnum(0,0,0)
                    if na==1:
                        motor[x1][i][0]+=d
                    if na==-1:
                        motor[x1][i][0]-=d
                if kina[x1][i]==0 and kinb[x1][i]==1:
                    nb=stepnum(0,0,0)
                    if nb==1:
                        motor[x1][i][1]-=d
                    if nb==-1:
                        motor[x1][i][1]+=d
            if kina[x1][i]==1 and kinb[x1][i]==1:
                if pan1!=0:
                    dx1=(motor[x1][i][0]-motor[x1][i][1])
                    for x2 in range(pair):
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(nu[x2]):
                            if kina[x2][ton1]==1 and kinb[x2][ton1]==1:
                                dxc[x2][ton1]=motor[x2][ton1][0]-motor[x2][ton1][1]
                                if x1==x2 and i==ton1:
                                    dxc[x2][ton1]+=d
                                if dxc[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                if dxc[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                            if kina[x2][ton1]==0 or kinb[x2][ton1]==0:
                                dxc[x2][ton1]=0
                        for ton1 in range(nu1[x2]):
                            if kinap[x2][ton1]==1 and kinbp[x2][ton1]==1:
                                dxcp[x2][ton1]=motorp[x2][ton1][0]-motorp[x2][ton1][1]
                                if dxcp[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                if dxcp[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                            if kinap[x2][ton1]==0 or kinbp[x2][ton1]==0:
                                dxcp[x2][ton1]=0
                        for ton1 in range(nu1x[x2]):
                            if kinapx[x2][ton1]==1 and kinbpx[x2][ton1]==1:
                                dxcpx[x2][ton1]=motorpx[x2][ton1][0]-motorpx[x2][ton1][1]
                                if dxcpx[x2][ton1]>xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                if dxcpx[x2][ton1]<-xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                            if kinapx[x2][ton1]==0 or kinbpx[x2][ton1]==0:
                                dxcpx[x2][ton1]=0
                    for ton1 in range(Nma[0]):
                        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                            F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                            F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                    for ton1 in range(Nma[1]):
                        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                            F3[0]-=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                            F3[1]+=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                    for ton1 in range(Nma1[0]):
                        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                            F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                            F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                    for ton1 in range(Nma1[1]):
                        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                            F4[0]-=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                            F4[1]+=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                                                #增加为负，减少为正
                    zj0[x1][i],zj1[x1][i],zj2[x1][i],zj3[x1][i],zj4[x1][i],zj5[x1][i],zj6[x1][i],zj7[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxc[x1][i]=dx1-d
                    if dx1+d>xdmotor:
                        F3[x1]-=(dx1+d-xdmotor)*Km
                        F2[x1]-=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F3[x1]-=(dx1+d+xdmotor)*Km
                        F2[x1]-=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F3[x1]+=(dx1-d-xdmotor)*Km
                        F2[x1]+=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F3[x1]+=(dx1-d+xdmotor)*Km
                        F2[x1]+=(dx1-d+xdmotor)*Km
                    js0[x1][i],js1[x1][i],js2[x1][i],js3[x1][i],js4[x1][i],js5[x1][i],js6[x1][i],js7[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl1[x1][i],ypsl2[x1][i]=energy(x1,i)

                na=stepnum(ypsl0,ypsl1[x1][i],ypsl2[x1][i])
                nb=stepnum(ypsl0,ypsl1[x1][i],ypsl2[x1][i])
                if motor[x1][i][0]+d>l2r[x1] and na==1:
                    na=0
                if motor[x1][i][1]-d<l3l[x1] and nb==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na>0.9:
                    zq=0
                    for xn in range(nu[x1]):
                        if motor[x1][i][0]+d>=motor[x1][xn][0]-0.1 and  motor[x1][i][0]<motor[x1][xn][0] and kina[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motor[x1][i][0]+=d
                        pan1=1
                        pan[x1][i]=1
                if na<-0.9:
                    zq=0
                    for xn in range(nu[x1]):
                        if motor[x1][i][0]-d<=motor[x1][xn][0] and motor[x1][i][0]>motor[x1][xn][0] and kina[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motor[x1][i][0]-=d
                        pan1=1
                        pan[x1][i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu[x1]):
                        if motor[x1][i][1]-d<=motor[x1][xn][1] and  motor[x1][i][1]>motor[x1][xn][1] and kinb[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motor[x1][i][1]-=d
                        pan1=1
                        pan[x1][i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu[x1]):
                        if motor[x1][i][1]+d>=motor[x1][xn][1] and  motor[x1][i][1]<motor[x1][xn][1] and kinb[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motor[x1][i][1]+=d
                        pan1=1
                        pan[x1][i]=1
                if pan[x1][i]==1 and zq==0:
                    #print('mz',ypsl0,ypsl1[i],ypsl2[i])
                    if na>0.1 or nb>0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=zj1[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=zj2[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=zj3[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=zj2[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=zj4[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=zj1[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=zj3[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=zj1[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=zj2[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=zj4[x1][i][x2]
                            l2l[x2]+=zj1[x1][i][x2]
                            l2r[x2]+=zj1[x1][i][x2]
                            l3l[x2]+=zj2[x1][i][x2]
                            l3r[x2]+=zj2[x1][i][x2]
                            l1l[x2]+=zj3[x1][i][x2]
                            l1r[x2]+=zj3[x1][i][x2]
                            l4l[x2]+=zj4[x1][i][x2]
                            l4r[x2]+=zj4[x1][i][x2]
                            Ximt2[x2]+=zj1[x1][i][x2]
                            Ximt3[x2]+=zj2[x1][i][x2]
                            Xkmt1[x2]+=zj3[x1][i][x2]
                            Xkmt4[x2]+=zj4[x1][i][x2]
                        Xpole+=zj0[x1][i]
                        Xpo+=zj7[x1][i]
                        Xkin+=zj5[x1][i]
                        Xkin1+=zj6[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=zj1[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=zj1[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=zj2[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=zj2[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=zj3[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=zj3[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=zj4[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=zj4[x1][i][1]
                    if na<-0.1 or nb<-0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=js1[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=js2[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=js3[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=js2[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=js4[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=js1[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=js3[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=js1[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=js2[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=js4[x1][i][x2]
                            l2l[x2]+=js1[x1][i][x2]
                            l2r[x2]+=js1[x1][i][x2]
                            l3l[x2]+=js2[x1][i][x2]
                            l3r[x2]+=js2[x1][i][x2]
                            l1l[x2]+=js3[x1][i][x2]
                            l1r[x2]+=js3[x1][i][x2]
                            l4l[x2]+=js4[x1][i][x2]
                            l4r[x2]+=js4[x1][i][x2]
                            Ximt2[x2]+=js1[x1][i][x2]
                            Ximt3[x2]+=js2[x1][i][x2]
                            Xkmt1[x2]+=js3[x1][i][x2]
                            Xkmt4[x2]+=js4[x1][i][x2]
                        Xpole+=js0[x1][i]
                        Xpo+=js7[x1][i]
                        Xkin+=js5[x1][i]
                        Xkin1+=js6[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=js1[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=js1[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=js2[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=js2[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=js3[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=js3[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=js4[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=js4[x1][i][1]
                    ypsl0=energy0()
                    ypsl0p=ypsl0
                    ypsl0px=ypsl0
                    ckF,de10,de11,de12,de13,de14,de15,de16,de17=checkF()
                    if ckF==1:
                        print(kk)
                        
                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        pan[x1][i]=0


    for x1 in range(pair):
        for i in range(nu1[x1]):
            if kinap[x1][i]+kinbp[x1][i]<2:
                if pan1==1:
                    panp[x1][i]=0
                if kinap[x1][i]==1 and kinbp[x1][i]==0:
                    na=stepnum(0,0,0)
                    if na==1:
                        motorp[x1][i][0]+=d
                    if na==-1:
                        motorp[x1][i][0]-=d
                if kinap[x1][i]==0 and kinbp[x1][i]==1:
                    nb=stepnum(0,0,0)
                    if nb==1:
                        motorp[x1][i][1]-=d
                    if nb==-1:
                        motorp[x1][i][1]+=d
            if kinap[x1][i]==1 and kinbp[x1][i]==1:
                if pan1!=0:
                    dx1=(motorp[x1][i][0]-motorp[x1][i][1])
                    for x2 in range(pair):
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(nu[x2]):
                            if kina[x2][ton1]==1 and kinb[x2][ton1]==1:
                                dxc[x2][ton1]=motor[x2][ton1][0]-motor[x2][ton1][1]

                                if dxc[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                if dxc[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                            if kina[x2][ton1]==0 or kinb[x2][ton1]==0:
                                dxc[x2][ton1]=0
                        for ton1 in range(nu1[x2]):
                            if kinap[x2][ton1]==1 and kinbp[x2][ton1]==1:
                                dxcp[x2][ton1]=motorp[x2][ton1][0]-motorp[x2][ton1][1]
                                if x1==x2 and i==ton1:
                                    dxcp[x2][ton1]+=d
                                if dxcp[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                if dxcp[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                            if kinap[x2][ton1]==0 or kinbp[x2][ton1]==0:
                                dxcp[x2][ton1]=0
                        for ton1 in range(nu1x[x2]):
                            if kinapx[x2][ton1]==1 and kinbpx[x2][ton1]==1:
                                dxcpx[x2][ton1]=motorpx[x2][ton1][0]-motorpx[x2][ton1][1]
                                if dxcpx[x2][ton1]>xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                if dxcpx[x2][ton1]<-xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                            if kinapx[x2][ton1]==0 or kinbpx[x2][ton1]==0:
                                dxcpx[x2][ton1]=0
                    for ton1 in range(Nma[0]):
                        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                            F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                            F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                    for ton1 in range(Nma[1]):
                        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                            F3[0]-=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                            F3[1]+=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                    for ton1 in range(Nma1[0]):
                        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                            F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                            F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                    for ton1 in range(Nma1[1]):
                        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                            F4[0]-=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                            F4[1]+=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                    zj0p[x1][i],zj1p[x1][i],zj2p[x1][i],zj3p[x1][i],zj4p[x1][i],zj5p[x1][i],zj6p[x1][i],zj7p[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxcp[x1][i]=dx1-d
                    if dx1+d>xdmotor:
                        F3[x1]-=(dx1+d-xdmotor)*Km
                        F1[x1]-=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F3[x1]-=(dx1+d+xdmotor)*Km
                        F1[x1]-=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F3[x1]+=(dx1-d-xdmotor)*Km
                        F1[x1]+=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F3[x1]+=(dx1-d+xdmotor)*Km
                        F1[x1]+=(dx1-d+xdmotor)*Km
                    
                    js0p[x1][i],js1p[x1][i],js2p[x1][i],js3p[x1][i],js4p[x1][i],js5p[x1][i],js6p[x1][i],js7p[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl1p[x1][i],ypsl2p[x1][i]=energyp(x1,i)

                na=stepnum(ypsl0p,ypsl1p[x1][i],ypsl2p[x1][i])
                nb=stepnum(ypsl0p,ypsl1p[x1][i],ypsl2p[x1][i])

                if motorp[x1][i][0]+d>l1r[x1] and na==1:
                    na=0
                if motorp[x1][i][1]-d<l3l[x1] and nb==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na>0.9:
                    zq=0
                    for xn in range(nu1[x1]):
                        if motorp[x1][i][0]+d>=motorp[x1][xn][0] and  motorp[x1][i][0]<motorp[x1][xn][0] and kinap[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorp[x1][i][0]+=d
                        pan1=1
                        panp[x1][i]=1
                if na<-0.9:
                    zq=0
                    for xn in range(nu1[x1]):
                        if motorp[x1][i][0]-d<=motorp[x1][xn][0] and motorp[x1][i][0]>motorp[x1][xn][0] and kinap[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorp[x1][i][0]-=d
                        pan1=1
                        panp[x1][i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu1[x1]):
                        if motorp[x1][i][1]-d<=motorp[x1][xn][1] and  motorp[x1][i][1]>motorp[x1][xn][1] and kinbp[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorp[x1][i][1]-=d
                        pan1=1
                        panp[x1][i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu1[x1]):
                        if motorp[x1][i][1]+d>=motorp[x1][xn][1] and  motorp[x1][i][1]<motorp[x1][xn][1] and kinbp[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorp[x1][i][1]+=d
                        pan1=1
                        panp[x1][i]=1
                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        panp[x1][i]=0
                if panp[x1][i]==1 and zq==0:
                    #print('mpz',ypsl0p,ypsl1p[i],ypsl2p[i],na,nb,motorp[i][0]-motorp[i][1])
                    if na>0.1 or nb>0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=zj1p[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=zj2p[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=zj3p[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=zj2p[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=zj4p[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=zj1p[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=zj3p[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=zj1p[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=zj2p[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=zj4p[x1][i][x2]
                            l2l[x2]+=zj1p[x1][i][x2]
                            l2r[x2]+=zj1p[x1][i][x2]
                            l3l[x2]+=zj2p[x1][i][x2]
                            l3r[x2]+=zj2p[x1][i][x2]
                            l1l[x2]+=zj3p[x1][i][x2]
                            l1r[x2]+=zj3p[x1][i][x2]
                            l4l[x2]+=zj4p[x1][i][x2]
                            l4r[x2]+=zj4p[x1][i][x2]
                            Ximt2[x2]+=zj1p[x1][i][x2]
                            Ximt3[x2]+=zj2p[x1][i][x2]
                            Xkmt1[x2]+=zj3p[x1][i][x2]
                            Xkmt4[x2]+=zj4p[x1][i][x2]
                        Xpole+=zj0p[x1][i]
                        Xpo+=zj7p[x1][i]
                        Xkin+=zj5p[x1][i]
                        Xkin1+=zj6p[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=zj1p[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=zj1p[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=zj2p[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=zj2p[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=zj3p[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=zj3p[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=zj4p[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=zj4p[x1][i][1]
                    if na<-0.1 or nb<-0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=js1p[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=js2p[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=js3p[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=js2p[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=js4p[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=js1p[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=js3p[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=js1p[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=js2p[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=js4p[x1][i][x2]
                            l2l[x2]+=js1p[x1][i][x2]
                            l2r[x2]+=js1p[x1][i][x2]
                            l3l[x2]+=js2p[x1][i][x2]
                            l3r[x2]+=js2p[x1][i][x2]
                            l1l[x2]+=js3p[x1][i][x2]
                            l1r[x2]+=js3p[x1][i][x2]
                            l4l[x2]+=js4p[x1][i][x2]
                            l4r[x2]+=js4p[x1][i][x2]
                            Ximt2[x2]+=js1p[x1][i][x2]
                            Ximt3[x2]+=js2p[x1][i][x2]
                            Xkmt1[x2]+=js3p[x1][i][x2]
                            Xkmt4[x2]+=js4p[x1][i][x2]
                        Xpole+=js0p[x1][i]
                        Xpo+=js7p[x1][i]
                        Xkin+=js5p[x1][i]
                        Xkin1+=js6p[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=js1p[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=js1p[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=js2p[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=js2p[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=js3p[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=js3p[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=js4p[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=js4p[x1][i][1]
                    ypsl0=energy0()
                    ypsl0p=ypsl0
                    ypsl0px=ypsl0
                    

    for x1 in range(pair):
        for i in range(nu1x[x1]):
            if kinapx[x1][i]+kinbpx[x1][i]<2:
                panpx[x1][i]=0
                if kinapx[x1][i]==1 and kinbpx[x1][i]==0:
                    na=stepnum(0,0,0)
                    if na==1:
                        motorpx[x1][i][0]-=d
                    if na==-1:
                        motorpx[x1][i][0]+=d
                if kinapx[x1][i]==0 and kinbpx[x1][i]==1:
                    nb=stepnum(0,0,0)
                    if nb==1:
                        motorpx[x1][i][1]+=d
                    if nb==-1:
                        motorpx[x1][i][1]-=d
            if kinapx[x1][i]==1 and kinbpx[x1][i]==1:
                if pan1!=0:
                    dx1=(motorpx[x1][i][0]-motorpx[x1][i][1])
                    for x2 in range(pair):
                        F3[x2]=Kp2*(Xpo-l3r[x2])
                        F1[x2]=Kp1*(l1l[x2]-Xpole)+Kp3*(l1r[x2]-Xkin)
                        F2[x2]=Kp2*(l2l[x2]-Xpole)
                        F4[x2]=Kp1*(Xpo-l4r[x2])+Kp3*(Xkin1-l4l[x2])
                        for ton1 in range(nu[x2]):
                            if kina[x2][ton1]==1 and kinb[x2][ton1]==1:
                                dxc[x2][ton1]=motor[x2][ton1][0]-motor[x2][ton1][1]
                                if dxc[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]-xdmotor)*Km
                                if dxc[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                                    F2[x2]+=(dxc[x2][ton1]+xdmotor)*Km
                            if kina[x2][ton1]==0 or kinb[x2][ton1]==0:
                                dxc[x2][ton1]=0
                        for ton1 in range(nu1[x2]):
                            if kinap[x2][ton1]==1 and kinbp[x2][ton1]==1:
                                dxcp[x2][ton1]=motorp[x2][ton1][0]-motorp[x2][ton1][1]

                                if dxcp[x2][ton1]>xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]-xdmotor)*Km
                                if dxcp[x2][ton1]<-xdmotor:
                                    F3[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                                    F1[x2]+=(dxcp[x2][ton1]+xdmotor)*Km
                            if kinap[x2][ton1]==0 or kinbp[x2][ton1]==0:
                                dxcp[x2][ton1]=0
                        for ton1 in range(nu1x[x2]):
                            if kinapx[x2][ton1]==1 and kinbpx[x2][ton1]==1:
                                dxcpx[x2][ton1]=motorpx[x2][ton1][0]-motorpx[x2][ton1][1]
                                if x1==x2 and i==ton1:
                                    dxcpx[x2][ton1]+=d
                                if dxcpx[x2][ton1]>xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]-xdmotor)*Km
                                if dxcpx[x2][ton1]<-xdmotor:
                                    F4[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                                    F2[x2]-=(dxcpx[x2][ton1]+xdmotor)*Km
                            if kinapx[x2][ton1]==0 or kinbpx[x2][ton1]==0:
                                dxcpx[x2][ton1]=0
                    for ton1 in range(Nma[0]):
                        if kinam[0][ton1]==1 and kinbm[0][ton1]==1:
                            F2[0]+=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                            F2[1]-=Kpn5*(numa[0][ton1][0]-numa[0][ton1][1])
                    for ton1 in range(Nma[1]):
                        if kinam[1][ton1]==1 and kinbm[1][ton1]==1:
                            F3[0]-=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                            F3[1]+=Kpn5*(numa[1][ton1][0]-numa[1][ton1][1])
                    for ton1 in range(Nma1[0]):
                        if kinam1[0][ton1]==1 and kinbm1[0][ton1]==1:
                            F1[0]+=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                            F1[1]-=Kpn5*(numa1[0][ton1][0]-numa1[0][ton1][1])
                    for ton1 in range(Nma1[1]):
                        if kinam1[1][ton1]==1 and kinbm1[1][ton1]==1:
                            F4[0]-=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                            F4[1]+=Kpn5*(numa1[1][ton1][0]-numa1[1][ton1][1])
                    zj0px[x1][i],zj1px[x1][i],zj2px[x1][i],zj3px[x1][i],zj4px[x1][i],zj5px[x1][i],zj6px[x1][i],zj7px[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    dxcpx[x1][i]=dx1-d
                    if dx1+d>xdmotor:
                        F4[x1]+=(dx1+d-xdmotor)*Km
                        F2[x1]+=(dx1+d-xdmotor)*Km
                    if dx1+d<-xdmotor:
                        F4[x1]+=(dx1+d+xdmotor)*Km
                        F2[x1]+=(dx1+d+xdmotor)*Km
                    if dx1-d>xdmotor:
                        F4[x1]-=(dx1-d-xdmotor)*Km
                        F2[x1]-=(dx1-d-xdmotor)*Km
                    if dx1-d<-xdmotor:
                        F4[x1]-=(dx1-d+xdmotor)*Km
                        F2[x1]-=(dx1-d+xdmotor)*Km
                    js0px[x1][i],js1px[x1][i],js2px[x1][i],js3px[x1][i],js4px[x1][i],js5px[x1][i],js6px[x1][i],js7px[x1][i]=jiao(dxc,dxcp,dxcpx,F1,F2,F3,F4,monum,monump,monumpx)
                    ypsl1px[x1][i],ypsl2px[x1][i]=energypx(x1,i)

                na=stepnum(ypsl0px,ypsl2px[x1][i],ypsl1px[x1][i])
                nb=stepnum(ypsl0px,ypsl2px[x1][i],ypsl1px[x1][i])
                if motorpx[x1][i][0]-d<l4l[x1] and na==1:
                    na=0
                if motorpx[x1][i][1]+d>l2r[x1] and nb==1:
                    nb=0
                if abs(na)>0.9 and abs(nb)>0.9:
                    if na+nb>=-1:
                        na=min(1,round(na+nb))
                    if na+nb<-1:
                        na=-1
                    nb=0
                if na<-0.9:
                    zq=0
                    for xn in range(nu1x[x1]):
                        if motorpx[x1][i][0]+d>=motorpx[x1][xn][0] and  motorpx[x1][i][0]<motorpx[x1][xn][0] and kinapx[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[x1][i][0]+=d
                        pan1=1
                        panpx[x1][i]=1

                if na>0.9:
                    zq=0
                    for xn in range(nu1x[x1]):
                        if motorpx[x1][i][0]-d<=motorpx[x1][xn][0] and motorpx[x1][i][0]>motorpx[x1][xn][0] and kinapx[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[x1][i][0]-=d
                        pan1=1
                        panpx[x1][i]=1
                if nb<-0.9 and na==0:
                    zq=0
                    for xn in range(nu1x[x1]):
                        if motorpx[x1][i][1]-d<=motorpx[x1][xn][1] and  motorpx[x1][i][1]>motorpx[x1][xn][1] and kinbpx[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[x1][i][1]-=d
                        pan1=1
                        panpx[x1][i]=1
                if nb>0.9 and na==0:
                    zq=0
                    for xn in range(nu1x[x1]):
                        if motorpx[x1][i][1]+d>=motorpx[x1][xn][1] and  motorpx[x1][i][1]<motorpx[x1][xn][1] and kinbpx[x1][xn]==1:
                            zq=3
                    if zq==0:
                        motorpx[x1][i][1]+=d
                        pan1=1
                        panpx[x1][i]=1
                if panpx[x1][i]==1 and zq==0:
                    #print('mpxz',ypsl0px,ypsl2px[i],ypsl1px[i],dxcpx[i])
                    if na<-0.1 or nb<-0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=zj1px[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=zj2px[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=zj3px[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=zj2px[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=zj4px[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=zj1px[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=zj3px[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=zj1px[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=zj2px[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=zj4px[x1][i][x2]
                            l2l[x2]+=zj1px[x1][i][x2]
                            l2r[x2]+=zj1px[x1][i][x2]
                            l3l[x2]+=zj2px[x1][i][x2]
                            l3r[x2]+=zj2px[x1][i][x2]
                            l1l[x2]+=zj3px[x1][i][x2]
                            l1r[x2]+=zj3px[x1][i][x2]
                            l4l[x2]+=zj4px[x1][i][x2]
                            l4r[x2]+=zj4px[x1][i][x2]
                            Ximt2[x2]+=zj1px[x1][i][x2]
                            Ximt3[x2]+=zj2px[x1][i][x2]
                            Xkmt1[x2]+=zj3px[x1][i][x2]
                            Xkmt4[x2]+=zj4px[x1][i][x2]
                        Xpole+=zj0px[x1][i]
                        Xpo+=zj7px[x1][i]
                        Xkin+=zj5px[x1][i]
                        Xkin1+=zj6px[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=zj1px[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=zj1px[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=zj2px[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=zj2px[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=zj3px[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=zj3px[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=zj4px[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=zj4px[x1][i][1]
                    if na>0.1 or nb>0.1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=js1px[x1][i][x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=js2px[x1][i][x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=js3px[x1][i][x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=js2px[x1][i][x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=js4px[x1][i][x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=js1px[x1][i][x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=js3px[x1][i][x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=js1px[x1][i][x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=js2px[x1][i][x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=js4px[x1][i][x2]
                            l2l[x2]+=js1px[x1][i][x2]
                            l2r[x2]+=js1px[x1][i][x2]
                            l3l[x2]+=js2px[x1][i][x2]
                            l3r[x2]+=js2px[x1][i][x2]
                            l1l[x2]+=js3px[x1][i][x2]
                            l1r[x2]+=js3px[x1][i][x2]
                            l4l[x2]+=js4px[x1][i][x2]
                            l4r[x2]+=js4px[x1][i][x2]
                            Ximt2[x2]+=js1px[x1][i][x2]
                            Ximt3[x2]+=js2px[x1][i][x2]
                            Xkmt1[x2]+=js3px[x1][i][x2]
                            Xkmt4[x2]+=js4px[x1][i][x2]
                        Xpole+=js0px[x1][i]
                        Xpo+=js7px[x1][i]
                        Xkin+=js5px[x1][i]
                        Xkin1+=js6px[x1][i]
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=js1px[x1][i][0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=js1px[x1][i][1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=js2px[x1][i][0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=js2px[x1][i][1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=js3px[x1][i][0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=js3px[x1][i][1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=js4px[x1][i][0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=js4px[x1][i][1]
                    ypsl0=energy0()
                    ypsl0p=ypsl0
                    ypsl0px=ypsl0
                    ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                    if ckF==1:
                        print(kk)
                if pan1==1:
                    if abs(na)+abs(nb)==0:
                        panpx[x1][i]=0  

    if pan1==1:
        for x1 in range(pair):
            if sum(pan[x1])+sum(panp[x1])+sum(panpx[x1])!=0:
                pan1=1
                break
            pan1=0
#############################################             
############################################################
    for x1 in range(pair):
        tuogeng=0
        for xbb in range(nu[x1]):
            xpa=0
            if kina[x1][xbb]==0 and kinb[x1][xbb]==1:
                xpa=1
                kinb[x1][xbb]=dissociation(0,0,0,0)
                if motor[x1][xbb][1]>=l2l[x1] and motor[x1][xbb][1]<=l2r[x1]:
                    kina[x1][xbb]=rebind()
                if kinb[x1][xbb]==1 and kina[x1][xbb]==1:
                    lian1=0
                    for xn in range(nu[x1]):
                        if kina[x1][xn]==1 and xn!=xbb and abs(l2r[x1]-round((l2r[x1]-motor[x1][xbb][1])/d)*d-motor[x1][xn][0])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l2r[x1]+round((motor[x1][xbb][1]-l2r[x1])/d)*d>=motor[x1][xbb][1]:
                            la=(round((motor[x1][xbb][1]-l2r[x1])/d)-1)*d
                        if l2r[x1]+round((motor[x1][xbb][1]-l2r[x1])/d)*d<motor[x1][xbb][1]:
                            la=(round((motor[x1][xbb][1]-l2r[x1])/d)+1)*d
                        for xn in range(nu[x1]):
                            if kina[x1][xn]==1 and xn!=xbb and abs(l2r[x1]+la-motor[x1][xn][0])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motor[x1][xbb][0]=l2r[x1]+round((motor[x1][xbb][1]-l2r[x1])/d)*d
                    if lian1==1:
                        motor[x1][xbb][0]=l2r[x1]+la
                    if lian1==2:
                        kina[x1][xbb]=0
                    pan1=1
                    monum[x1]+=1
            if xpa==0 and kina[x1][xbb]==1 and kinb[x1][xbb]==0:
                xpa=1
                kina[x1][xbb]=dissociation(0,0,0,0)
                if motor[x1][xbb][0]>=l3l[x1] and motor[x1][xbb][0]<=l3r[x1]:
                    kinb[x1][xbb]=rebind()
                if kinb[x1][xbb]==1 and kina[x1][xbb]==1:
                    lian1=0
                    for xn in range(nu[x1]):
                        if kinb[x1][xn]==1 and xn!=xbb and abs(l3l[x1]+round((motor[x1][xbb][0]-l3l[x1])/d)*d-motor[x1][xn][1])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l3l[x1]+round((motor[x1][xbb][0]-l3l[x1])/d)*d>=motor[x1][xbb][0]:
                            la=(round((motor[x1][xbb][0]-l3l[x1])/d)-1)*d
                        if l3l[x1]+round((motor[x1][xbb][0]-l3l[x1])/d)*d<motor[x1][xbb][0]:
                            la=(round((motor[x1][xbb][0]-l3l[x1])/d)+1)*d
                        for xn in range(nu[x1]):
                            if kinb[x1][xn]==1 and xn!=xbb and abs(l3l[x1]+la-motor[x1][xn][1])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motor[x1][xbb][1]=l3l[x1]+round((motor[x1][xbb][0]-l3l[x1])/d)*d
                    if lian1==1:
                        motor[x1][xbb][1]=l3l[x1]+la
                    if lian1==2:
                        kinb[x1][xbb]=0
                    pan1=1
                    monum[x1]+=1
            if xpa==0 and kina[x1][xbb]==1 and kinb[x1][xbb]==1:
                F=0
                if motor[x1][xbb][0]-motor[x1][xbb][1]>xdmotor:
                    F=Km*(motor[x1][xbb][0]-motor[x1][xbb][1]-xdmotor)
                if motor[x1][xbb][0]-motor[x1][xbb][1]<-xdmotor:
                    F=Km*(motor[x1][xbb][0]-motor[x1][xbb][1]+xdmotor)
                kina[x1][xbb]=dissociation(F,ypsl0px,ypsl1[x1][xbb],ypsl2[x1][xbb])
                kinb[x1][xbb]=dissociation(F,ypsl0px,ypsl1[x1][xbb],ypsl2[x1][xbb])
                if motor[x1][xbb][0]>l2r[x1] or motor[x1][xbb][0]<l2l[x1]:
                    kina[x1][xbb]=0
                if motor[x1][xbb][1]<l3l[x1] or motor[x1][xbb][1]>l3r[x1]:
                    kinb[x1][xbb]=0
                if kina[x1][xbb]==0 or kinb[x1][xbb]==0:
                    monum[x1]-=1
                    pan1=1
                    ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                    if ckF==1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=de1[x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=de2[x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=de3[x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=de2[x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=de4[x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=de1[x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=de3[x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=de1[x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=de2[x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=de4[x2]
                            l2l[x2]+=de1[x2]
                            l2r[x2]+=de1[x2]
                            l3l[x2]+=de2[x2]
                            l3r[x2]+=de2[x2]
                            l1l[x2]+=de3[x2]
                            l1r[x2]+=de3[x2]
                            l4l[x2]+=de4[x2]
                            l4r[x2]+=de4[x2]
                            Ximt2[x2]+=de1[x2]
                            Ximt3[x2]+=de2[x2]
                            Xkmt1[x2]+=de3[x2]
                            Xkmt4[x2]+=de4[x2]
                        Xpole+=de0
                        Xpo+=de7
                        Xkin+=de5
                        Xkin1+=de6
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=de1[0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=de1[1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=de2[0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=de2[1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=de3[0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=de3[1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=de4[0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=de4[1]
                        ypsl0=energy0()
                        ypsl0p=ypsl0
                        ypsl0px=ypsl0
                   
            if kina[x1][xbb]==0 and kinb[x1][xbb]==0:
                tuogeng=1
######################################################
        if tuogeng==1:
            for xn in range(nu[x1]):
                if kina[x1][xn]==0 and kinb[x1][xn]==0:
                    tuogeng=0
                    del motor[x1][xn]
                    del kina[x1][xn]
                    del kinb[x1][xn]
                    del dxc[x1][xn]
                    del ypsl1[x1][xn]
                    del ypsl2[x1][xn]
                    del zj7[x1][xn]
                    del js7[x1][xn]
                    del zj6[x1][xn]
                    del js6[x1][xn]
                    del zj5[x1][xn]
                    del js5[x1][xn]
                    del zj4[x1][xn]
                    del js4[x1][xn]
                    del zj3[x1][xn]
                    del js3[x1][xn]
                    del zj2[x1][xn]
                    del js2[x1][xn]
                    del zj1[x1][xn]
                    del js1[x1][xn]
                    del zj0[x1][xn]
                    del js0[x1][xn]
                    del pan[x1][xn]
                    nu[x1]-=1
                    break 
##################################################################
        tuogeng=0
        nu[x1]=len(kina[x1])
        for xbb in range(nu1[x1]):
            xpa=0
            if kinap[x1][xbb]==0 and kinbp[x1][xbb]==1:
                xpa=1
                kinbp[x1][xbb]=dissociation(0,0,0,0)
                if motorp[x1][xbb][1]>=l3l[x1] and motorp[x1][xbb][1]<=l1r[x1]:
                    kinap[x1][xbb]=rebind()
                if kinbp[x1][xbb]==1 and kinap[x1][xbb]==1:
                    lian1=0
                    for xn in range(nu1[x1]):
                        if kinap[x1][xn]==1 and xn!=xbb and abs(l1l[x1]-round((l1l[x1]-motorp[x1][xbb][1])/d)*d-motorp[x1][xn][0])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l1l[x1]+round((motorp[x1][xbb][1]-l1l[x1])/d)*d>=motorp[x1][xbb][1]:
                            la=(round((motorp[x1][xbb][1]-l1l[x1])/d)-1)*d
                        if l1l[x1]+round((motorp[x1][xbb][1]-l1l[x1])/d)*d<motorp[x1][xbb][1]:
                            la=(round((motorp[x1][xbb][1]-l1l[x1])/d)+1)*d
                        for xn in range(nu1[x1]):
                            if kinap[x1][xn]==1 and xn!=xbb and abs(l1l[x1]+la-motorp[x1][xn][0])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motorp[x1][xbb][0]=l1l[x1]+round((motorp[x1][xbb][1]-l1l[x1])/d)*d
                    if lian1==1:
                        motorp[x1][xbb][0]=l1l[x1]+la
                    if lian1==2:
                        kinap[x1][xbb]=0
                    pan1=1
                    monump[x1]+=1
            if xpa==0 and kinap[x1][xbb]==1 and kinbp[x1][xbb]==0:
                xpa=1
                kinap[x1][xbb]=dissociation(0,0,0,0)
                if motorp[x1][xbb][0]>=l3l[x1] and motorp[x1][xbb][0]<=l1r[x1]:
                    kinbp[x1][xbb]=rebind()
                if kinbp[x1][xbb]==1 and kinap[x1][xbb]==1:
                    lian1=0
                    for xn in range(nu1[x1]):
                        if kinbp[x1][xn]==1 and xn!=xbb and abs(l3l[x1]+round((motorp[x1][xbb][0]-l3l[x1])/d)*d-motorp[x1][xn][1])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l3l[x1]+round((motorp[x1][xbb][0]-l3l[x1])/d)*d>=motorp[x1][xbb][0]:
                            la=(round((motorp[x1][xbb][0]-l3l[x1])/d)-1)*d
                        if l3l[x1]+round((motorp[x1][xbb][0]-l3l[x1])/d)*d<motorp[x1][xbb][0]:
                            la=(round((motorp[x1][xbb][0]-l3l[x1])/d)+1)*d
                        for xn in range(nu1[x1]):
                            if kinbp[x1][xn]==1 and xn!=xbb and abs(l3l[x1]+la-motorp[x1][xn][1])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motorp[x1][xbb][1]=l3l[x1]+round((motorp[x1][xbb][0]-l3l[x1])/d)*d
                    if lian1==1:
                        motorp[x1][xbb][1]=l3l[x1]+la
                    if lian1==2:
                        kinbp[x1][xbb]=0
                    pan1=1
                    monump[x1]+=1
            if xpa==0 and kinap[x1][xbb]==1 and kinbp[x1][xbb]==1:
                F=0
                if motorp[x1][xbb][0]-motorp[x1][xbb][1]>xdmotor:
                    F=Km*(motorp[x1][xbb][0]-motorp[x1][xbb][1]-xdmotor)
                if motorp[x1][xbb][0]-motorp[x1][xbb][1]<-xdmotor:
                    F=Km*(motorp[x1][xbb][0]-motorp[x1][xbb][1]+xdmotor)
                kinap[x1][xbb]=dissociation(F,ypsl0px,ypsl1p[x1][xbb],ypsl2p[x1][xbb])
                kinbp[x1][xbb]=dissociation(F,ypsl0px,ypsl1p[x1][xbb],ypsl2p[x1][xbb])
                if motorp[x1][xbb][0]>l1r[x1] or motorp[x1][xbb][0]<l1l[x1]:
                    kinap[x1][xbb]=0
                if motorp[x1][xbb][1]<l3l[x1] or motorp[x1][xbb][1]>l3r[x1]:
                    kinbp[x1][xbb]=0
                if kinap[x1][xbb]==0 or kinbp[x1][xbb]==0:
                    monump[x1]-=1
                    pan1=1
                    ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                    if ckF==1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=de1[x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=de2[x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=de3[x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=de2[x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=de4[x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=de1[x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=de3[x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=de1[x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=de2[x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=de4[x2]
                            l2l[x2]+=de1[x2]
                            l2r[x2]+=de1[x2]
                            l3l[x2]+=de2[x2]
                            l3r[x2]+=de2[x2]
                            l1l[x2]+=de3[x2]
                            l1r[x2]+=de3[x2]
                            l4l[x2]+=de4[x2]
                            l4r[x2]+=de4[x2]
                            Ximt2[x2]+=de1[x2]
                            Ximt3[x2]+=de2[x2]
                            Xkmt1[x2]+=de3[x2]
                            Xkmt4[x2]+=de4[x2]
                        Xpole+=de0
                        Xpo+=de7
                        Xkin+=de5
                        Xkin1+=de6
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=de1[0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=de1[1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=de2[0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=de2[1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=de3[0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=de3[1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=de4[0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=de4[1]
                        ypsl0=energy0()
                        ypsl0p=ypsl0
                        ypsl0px=ypsl0
                    
            if kinap[x1][xbb]==0 and kinbp[x1][xbb]==0:
                tuogeng=1
        if tuogeng==1:
            for xn in range(nu1[x1]):
                if kinap[x1][xn]==0 and kinbp[x1][xn]==0:
                    tuogeng=0
                    del motorp[x1][xn]
                    del kinap[x1][xn]
                    del kinbp[x1][xn]
                    del dxcp[x1][xn]
                    del ypsl1p[x1][xn]
                    del ypsl2p[x1][xn]
                    del zj7p[x1][xn]
                    del js7p[x1][xn]
                    del zj6p[x1][xn]
                    del js6p[x1][xn]
                    del zj5p[x1][xn]
                    del js5p[x1][xn]
                    del zj4p[x1][xn]
                    del js4p[x1][xn]
                    del zj3p[x1][xn]
                    del js3p[x1][xn]
                    del zj2p[x1][xn]
                    del js2p[x1][xn]
                    del zj1p[x1][xn]
                    del js1p[x1][xn]
                    del zj0p[x1][xn]
                    del js0p[x1][xn]
                    del panp[x1][xn]
                    nu1[x1]-=1
                    break
        tuogeng=0
        for xbb in range(nu1x[x1]):
            xpa=0
            if kinapx[x1][xbb]==0 and kinbpx[x1][xbb]==1:
                xpa=1
                kinbpx[x1][xbb]=dissociation(0,0,0,0)
                if motorpx[x1][xbb][1]>=l4l[x1] and motorpx[x1][xbb][1]<=l2r[x1]:
                    kinapx[x1][xbb]=rebind()

                if kinbpx[x1][xbb]==1 and kinapx[x1][xbb]==1:
                    lian1=0
                    #print('111')
                    for xn in range(nu1x[x1]):
                        if kinapx[x1][xn]==1 and xn!=xbb and abs(l4l[x1]-round((l4l[x1]-motorpx[x1][xbb][1])/d)*d-motorpx[x1][xn][0])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l4l[x1]+round((motorpx[x1][xbb][1]-l4l[x1])/d)*d>=motorpx[x1][xbb][1]:
                            la=(round((motorpx[x1][xbb][1]-l4l[x1])/d)-1)*d
                        if l4l[x1]+round((motorpx[x1][xbb][1]-l4l[x1])/d)*d<motorpx[x1][xbb][1]:
                            la=(round((motorpx[x1][xbb][1]-l4l[x1])/d)+1)*d
                        for xn in range(nu1x[x1]):
                            if kinapx[x1][xn]==1 and xn!=xbb and abs(l4l[x1]+la-motorpx[x1][xn][0])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motorpx[x1][xbb][0]=l4l[x1]+round((motorpx[x1][xbb][1]-l4l[x1])/d)*d
                    if lian1==1:
                        motorpx[x1][xbb][0]=l4l[x1]+la
                    if lian1==2:
                        kinapx[x1][xbb]=0
                        #print('222')
                    pan1=1
                    monumpx[x1]+=1
            if xpa==0 and kinapx[x1][xbb]==1 and kinbpx[x1][xbb]==0:
                xpa=1
                kinapx[x1][xbb]=dissociation(0,0,0,0)
                if motorpx[x1][xbb][0]>=l4l[x1] and motorpx[x1][xbb][0]<=l2r[x1]:
                    kinbpx[x1][xbb]=rebind()
                if kinbpx[x1][xbb]==1 and kinapx[x1][xbb]==1:
                    lian1=0
                    for xn in range(nu1x[x1]):
                        if kinbpx[x1][xn]==1 and xn!=xbb and abs(l2l[x1]+round((motorpx[x1][xbb][0]-l2l[x1])/d)*d-motorpx[x1][xn][1])<1:
                            lian1=1
                            break
                    if lian1==1:
                        if l2l[x1]+round((motorpx[x1][xbb][0]-l2l[x1])/d)*d>=motorpx[x1][xbb][0]:
                            la=(round((motorpx[x1][xbb][0]-l2l[x1])/d)-1)*d
                        if l2l[x1]+round((motorpx[x1][xbb][0]-l2l[x1])/d)*d<motorpx[x1][xbb][0]:
                            la=(round((motorpx[x1][xbb][0]-l2l[x1])/d)+1)*d
                        for xn in range(nu1x[x1]):
                            if kinbpx[x1][xn]==1 and xn!=xbb and abs(l2l[x1]+la-motorpx[x1][xn][1])<1:
                                lian1=2
                                break
                    if lian1==0:
                        motorpx[x1][xbb][1]=l2l[x1]+round((motorpx[x1][xbb][0]-l2l[x1])/d)*d
                    if lian1==1:
                        motorpx[x1][xbb][1]=l2l[x1]+la
                    if lian1==2:
                        kinbpx[x1][xbb]=0
                    pan1=1
                    monumpx[x1]+=1
            if xpa==0 and kinapx[x1][xbb]==1 and kinbpx[x1][xbb]==1:
                F=0
                if motorpx[x1][xbb][0]-motorpx[x1][xbb][1]>xdmotor:
                    F=-Km*(motorpx[x1][xbb][0]-motorpx[x1][xbb][1]-xdmotor)
                if motorpx[x1][xbb][0]-motorpx[x1][xbb][1]<-xdmotor:
                    F=-Km*(motorpx[x1][xbb][0]-motorpx[x1][xbb][1]+xdmotor)
                kinapx[x1][xbb]=dissociation(F,ypsl0px,ypsl2px[x1][xbb],ypsl1px[x1][xbb])
                kinbpx[x1][xbb]=dissociation(F,ypsl0px,ypsl2px[x1][xbb],ypsl1px[x1][xbb])
                if motorpx[x1][xbb][0]>l4r[x1] or motorpx[x1][xbb][0]<l4l[x1]:
                    kinapx[x1][xbb]=0
                if motorpx[x1][xbb][1]<l2l[x1] or motorpx[x1][xbb][1]>l2r[x1]:
                    kinbpx[x1][xbb]=0
                if kinapx[x1][xbb]==0 or kinbpx[x1][xbb]==0:
                    monumpx[x1]-=1
                    pan1=1
                    ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                    if ckF==1:
                        for x2 in range(pair):
                            for xn in range(nu[x2]):
                                if kina[x2][xn]==1:
                                    motor[x2][xn][0]+=de1[x2]
                                if kinb[x2][xn]==1:
                                    motor[x2][xn][1]+=de2[x2]
                            for xn in range(nu1[x2]):
                                if kinap[x2][xn]==1:
                                    motorp[x2][xn][0]+=de3[x2]
                                if kinbp[x2][xn]==1:
                                    motorp[x2][xn][1]+=de2[x2]
                            for xn in range(nu1x[x2]):
                                if kinapx[x2][xn]==1:
                                    motorpx[x2][xn][0]+=de4[x2]
                                if kinbpx[x2][xn]==1:
                                    motorpx[x2][xn][1]+=de1[x2]
                            for xn in range(numc1[x2]):
                                motormc1[x2][xn]+=de3[x2]
                            for xn in range(numc2[x2]):
                                motormc2[x2][xn]+=de1[x2]
                            for xn in range(numc3[x2]):
                                motormc3[x2][xn]+=de2[x2]
                            for xn in range(numc4[x2]):
                                motormc4[x2][xn]+=de4[x2]
                            l2l[x2]+=de1[x2]
                            l2r[x2]+=de1[x2]
                            l3l[x2]+=de2[x2]
                            l3r[x2]+=de2[x2]
                            l1l[x2]+=de3[x2]
                            l1r[x2]+=de3[x2]
                            l4l[x2]+=de4[x2]
                            l4r[x2]+=de4[x2]
                            Ximt2[x2]+=de1[x2]
                            Ximt3[x2]+=de2[x2]
                            Xkmt1[x2]+=de3[x2]
                            Xkmt4[x2]+=de4[x2]
                        Xpole+=de0
                        Xpo+=de7
                        Xkin+=de5
                        Xkin1+=de6
                        for ton1 in range(Nma[0]):   
                            if kinam[0][ton1]==1:
                                numa[0][ton1][0]+=de1[0]
                            if kinbm[0][ton1]==1:
                                numa[0][ton1][1]+=de1[1]
                        for ton1 in range(Nma[1]):   
                            if kinam[1][ton1]==1:
                                numa[1][ton1][0]+=de2[0]
                            if kinbm[1][ton1]==1:
                                numa[1][ton1][1]+=de2[1]
                        for ton1 in range(Nma1[0]):   
                            if kinam1[0][ton1]==1:
                                numa1[0][ton1][0]+=de3[0]
                            if kinbm1[0][ton1]==1:
                                numa1[0][ton1][1]+=de3[1]
                        for ton1 in range(Nma1[1]):   
                            if kinam1[1][ton1]==1:
                                numa1[1][ton1][0]+=de4[0]
                            if kinbm1[1][ton1]==1:
                                numa1[1][ton1][1]+=de4[1]
                        ypsl0=energy0()
                        ypsl0p=ypsl0
                        ypsl0px=ypsl0
                        ckF,de0,de1,de2,de3,de4,de5,de6,de7=checkF()
                        if ckF==1:
                            for x2 in range(pair):
                                for xn in range(nu[x2]):
                                    if kina[x2][xn]==1:
                                        motor[x2][xn][0]+=de1[x2]
                                    if kinb[x2][xn]==1:
                                        motor[x2][xn][1]+=de2[x2]
                                for xn in range(nu1[x2]):
                                    if kinap[x2][xn]==1:
                                        motorp[x2][xn][0]+=de3[x2]
                                    if kinbp[x2][xn]==1:
                                        motorp[x2][xn][1]+=de2[x2]
                                for xn in range(nu1x[x2]):
                                    if kinapx[x2][xn]==1:
                                        motorpx[x2][xn][0]+=de4[x2]
                                    if kinbpx[x2][xn]==1:
                                        motorpx[x2][xn][1]+=de1[x2]
                                for xn in range(numc1[x2]):
                                    motormc1[x2][xn]+=de3[x2]
                                for xn in range(numc2[x2]):
                                    motormc2[x2][xn]+=de1[x2]
                                for xn in range(numc3[x2]):
                                    motormc3[x2][xn]+=de2[x2]
                                for xn in range(numc4[x2]):
                                    motormc4[x2][xn]+=de4[x2]
                                l2l[x2]+=de1[x2]
                                l2r[x2]+=de1[x2]
                                l3l[x2]+=de2[x2]
                                l3r[x2]+=de2[x2]
                                l1l[x2]+=de3[x2]
                                l1r[x2]+=de3[x2]
                                l4l[x2]+=de4[x2]
                                l4r[x2]+=de4[x2]
                                Ximt2[x2]+=de1[x2]
                                Ximt3[x2]+=de2[x2]
                                Xkmt1[x2]+=de3[x2]
                                Xkmt4[x2]+=de4[x2]
                            Xpole+=de0
                            Xpo+=de7
                            Xkin+=de5
                            Xkin1+=de6
                            for ton1 in range(Nma[0]):   
                                if kinam[0][ton1]==1:
                                    numa[0][ton1][0]+=de1[0]
                                if kinbm[0][ton1]==1:
                                    numa[0][ton1][1]+=de1[1]
                            for ton1 in range(Nma[1]):   
                                if kinam[1][ton1]==1:
                                    numa[1][ton1][0]+=de2[0]
                                if kinbm[1][ton1]==1:
                                    numa[1][ton1][1]+=de2[1]
                            for ton1 in range(Nma1[0]):   
                                if kinam1[0][ton1]==1:
                                    numa1[0][ton1][0]+=de3[0]
                                if kinbm1[0][ton1]==1:
                                    numa1[0][ton1][1]+=de3[1]
                            for ton1 in range(Nma1[1]):   
                                if kinam1[1][ton1]==1:
                                    numa1[1][ton1][0]+=de4[0]
                                if kinbm1[1][ton1]==1:
                                    numa1[1][ton1][1]+=de4[1]
                            ypsl0=energy0()
                            ypsl0p=ypsl0
                            ypsl0px=ypsl0
            if kinapx[x1][xbb]==0 and kinbpx[x1][xbb]==0:
                tuogeng=1
        if tuogeng==1:
            for xn in range(nu1x[x1]):
                if kinapx[x1][xn]==0 and kinbpx[x1][xn]==0:
                    tuogeng=0
                    del motorpx[x1][xn]
                    del kinapx[x1][xn]
                    del kinbpx[x1][xn]
                    del dxcpx[x1][xn]
                    del ypsl1px[x1][xn]
                    del ypsl2px[x1][xn]
                    del zj7px[x1][xn]
                    del js7px[x1][xn]
                    del zj6px[x1][xn]
                    del js6px[x1][xn]
                    del zj5px[x1][xn]
                    del js5px[x1][xn]
                    del zj4px[x1][xn]
                    del js4px[x1][xn]
                    del zj3px[x1][xn]
                    del js3px[x1][xn]
                    del zj2px[x1][xn]
                    del js2px[x1][xn]
                    del zj1px[x1][xn]
                    del js1px[x1][xn]
                    del zj0px[x1][xn]
                    del js0px[x1][xn]
                    del panpx[x1][xn]
                    nu1x[x1]=len(kinapx[x1])
                    break

    if  Xpole>0:
        break
    if a%1000==0:
        ads+=1
        print(round(a*h),round(numpy.mean(monum)),round(numpy.mean(monump)),round(numpy.mean(monumpx)),round(numpy.mean(nunum1)),round(Xkin,2),round(Xkin1,2),round(numpy.mean(lxian[0]),1),round(numpy.mean(lxian[1]),1),round(Xpo-Xpole),round(((ximta2[0])-Ximt2[0]),2),round(((ximta3[0])-Ximt3[0]),2),round(((ximta2[1])-Ximt2[1]),2),round(((xkmta1[0])-Xkmt1[0]),2),round(((xkmta4[0])-Xkmt4[0]),2),round(((xkmta1[1])-Xkmt1[1]),2),round(((xkmta4[1])-Xkmt4[1]),2),round((numpy.mean(desi2)-numpy.mean(desia2)),2),round((numpy.mean(desk1)-numpy.mean(deska1)),2),round((numpy.mean(pos1)-numpy.mean(posa1)),2),round((numpy.mean(pos2)-numpy.mean(posa2)),2),round(Kp2*(l2l[0]-Xpole),2),round(Kp2*(l2l[1]-Xpole),2))
        for x2 in range(pair):
            desia2[x2]=desi2[x2]
            deska1[x2]=desk1[x2]
            posa1[x2]=pos1[x2]
            posa2[x2]=pos2[x2]
            ximta2[x2]=Ximt2[x2]
            xkmta1[x2]=Xkmt1[x2]
            ximta3[x2]=Ximt3[x2]
            xkmta4[x2]=Xkmt4[x2]




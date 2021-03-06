# -*-coding:utf-8-*-
import numpy as np
import matplotlib.pyplot as plt

# 输入层数、折射率、厚度、入射角度
def inputParameter():            # 输入层数、折射率、厚度、入射角度
    k=int(input('请输入膜层层数：'))
    o=2*np.pi/360*float(input('请输入入射角度（0°≤θ≤90°）：'))
    i=0;n=[];d=[0]
    while i <k+2:
        if i == 0 :
            n0=input('请输入入射介质的折射率n0：')
            n.append(float(n0))
        elif i>0 and i<=k:
            print('请输入第',i,'层的折射率n',i,':',)
            nj=input()
            dj=float(input('厚度(nm)：'))
            n.append(float(nj))
            d.append(dj)
        else:
            ns = input('请输入基底介质的折射率ns：')
            n.append(float(ns))
            d.append(0)
        i=i+1
    return n,d,o,k
# 要计算波长范围输入
def wavelengthInput():            # 要计算波长范围输入
    wavelengthMin=int(input('计算波长的最少值'))
    spacing=float(input('计算波长的间隔'))
    wavelengthMax=int(input('计算波长的最大值'))
    Wavelength=np.arange(wavelengthMin,wavelengthMax,spacing)     #产生波长序列
    return Wavelength
# cosφk的计算
def cos_RefractionAngle_k(o,n):        # cosφk的计算
    cos_RA=[]
    sin_RAe_0=np.sin(o)                      # sinφ0 的计算
    for i in n:
        sin_RAe_k=n[0]*sin_RAe_0/i
        cos_RefractionAngle=(1-sin_RAe_k**2)**0.5             #  cosφk=(1-(n0*sinφ0 / nk)²)^0.5
        cos_RA.append(cos_RefractionAngle)
    return cos_RA
#相位厚度
def PhaseThickness_k(d,cos_ar_k,wavelength,k,n):      #相位厚度
    ptk=[];kk=0
    while kk<=k+1:
        if kk==0:
            ptk.append(0)                                             #入射介质相位厚度
        else:
            PT_k=2*np.pi*cos_ar_k[kk]*n[kk]*d[kk]/wavelength      #2π*dk*cosφk / λ
            ptk.append(PT_k)
        kk=kk+1
    return ptk
# p分量的光学导纳
def Pyk(n,cos_ar_k,k):     # p分量的光学导纳
    hdsj=0;pyk=[]
    while hdsj<=k+1:
        pyk.append(n[hdsj]/cos_ar_k[hdsj])                   #  ηk=nk / cosφk
        hdsj = hdsj + 1
    return pyk
# s分量的光学导纳
def Syk(n,cos_ar_k,k):     # s分量的光学导纳
    hdsj1 = 0;syk = []
    while hdsj1 <= k + 1:
        syk.append(n[hdsj1] * cos_ar_k[hdsj1])              #  ηk=nk*cosφk
        hdsj1=hdsj1+1
    return syk
#菲涅尔系数r_s
def FresnelCoefficient_r_s(cos_ar_k,n,k):
    r_s = []
    for i in range(k + 1):
        r_s1=(n[i]*cos_ar_k[i])-(n[i+1]*cos_ar_k[i+1])
        r_s2=(n[i]*cos_ar_k[i])+(n[i+1]*cos_ar_k[i+1])
        r_s3=r_s1/r_s2
        r_s.append(r_s3)
    return r_s
#菲涅尔系数r_p
def FresnelCoefficient_r_p(cos_ar_k,n,k):
    r_p = []
    for i in range(k + 1):
        r_p1 = (n[i] * cos_ar_k[i+1]) - (n[i + 1] * cos_ar_k[i])
        r_p2 = (n[i] * cos_ar_k[i+ 1]) + (n[i + 1] * cos_ar_k[i])
        r_p3 = r_p1 / r_p2
        r_p.append(r_p3)
    return r_p
# 递推法
def RecurrenceMethod(xwhd,r_ps,k):
    r_ps = r_ps[::-1]  # 反转列表
    xwhd = xwhd[::-1]
    hjsahj = []
    hjsahj_4 = 0
    for i in range(k +1):
        if i == 1:
            hjsahj_1 = r_ps[i] + r_ps[i - 1] * np.exp(-1j * 2 * xwhd[i])
            hjsahj_2 = 1 + r_ps[i] * r_ps[i - 1] * np.exp(-1j * 2 * xwhd[i])
            hjsahj_3 = hjsahj_1 / hjsahj_2
            hjsahj_4 = hjsahj_3
            hjsahj.append(hjsahj_3)
        elif i > 1:
            hjsahj_1 = r_ps[i] + hjsahj_4 * np.exp(-1j * 2 * xwhd[i])
            hjsahj_2 = 1 + r_ps[i] * hjsahj_4 * np.exp(-1j * 2 * xwhd[i])
            hjsahj_3 = hjsahj_1 / hjsahj_2
            hjsahj_4 = hjsahj_3
            hjsahj.append(hjsahj_3)
    jhgd = hjsahj[-1] * (hjsahj[-1].conjugate())
    return jhgd.real




n,d,o,k=inputParameter()            # 输入层数k、折射率n、厚度d、入射角度o
Wavelength=wavelengthInput()     #要计算波长范围输入
reflectivity=[]           #能量反射率列表，空
for wavelength in Wavelength:
    cos_ar_k=cos_RefractionAngle_k(o,n)      #  cosφk的计算
    xwhd=PhaseThickness_k(d,cos_ar_k,wavelength,k,n)   #相位厚度
    PYK=Pyk(n,cos_ar_k,k)                     # p分量的光学导纳
    SYK=Syk(n,cos_ar_k,k)                    # s分量的光学导纳
    r_p=FresnelCoefficient_r_p(cos_ar_k,n,k)
    r_s=FresnelCoefficient_r_s(cos_ar_k,n,k)
    ppp=RecurrenceMethod(xwhd,r_p,k)
    sss=RecurrenceMethod(xwhd,r_s,k)
    rrrr=(ppp+sss)/2
    reflectivity.append(rrrr * 100)
print(reflectivity)
plt.plot(Wavelength,reflectivity,'r')
plt.xlabel('Wavelength (nm)')      # 设置 x轴的名称
plt.ylabel('Reflectivity (%)')               # 设置 y轴的名称
plt.xlim(380,780)                # 设置x轴的坐标的范围
#plt.ylim(0,5)                  # 设置y轴的坐标的范围
plt.show()

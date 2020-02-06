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
# 矩阵[BC]计算
def CharacteristicMatrix(YK,xwhd,k):    #矩阵[BC]
    bc_2=[1,YK[-1]]
    gdsh=0;bc_1=[]
    while gdsh<=k:                       #各层特征矩阵连乘
        if gdsh==0:
            a0=np.identity(2)               #创建单位矩阵
            bc_1.append(a0)
        else:
            m11=np.cos(xwhd[gdsh]);m22=np.cos(xwhd[gdsh])          #m11,m12,m21,m22特征矩阵的各元素
            m12=1j*np.sin(xwhd[gdsh])/YK[gdsh]
            m21=1j*np.sin(xwhd[gdsh])*YK[gdsh]
            a1=np.array([[m11,m12],[m21,m22]])                 #创建特征矩阵
            a2=np.dot(bc_1[-1],a1)                                       #前一个特征矩阵与后一个特征矩阵
            bc_1.append(a2)
        gdsh=gdsh+1
    bc_11=bc_1[-1]
    bc=np.dot(bc_11,bc_2)                                                            #得到BC矩阵
    return bc
# 衬底和膜系组合的能量反射率
def reflex(BCp,BCs,PYK_0,SYK_0):                     # 衬底和膜系组合的能量反射率
    py=BCp[1]/BCp[0]           # p的等效导纳
    sy=BCs[1]/BCs[0]            # s的等效导纳
    prj=(PYK_0-py)/(PYK_0+py)        #p振幅反射率
    srj=(SYK_0-sy)/(SYK_0+sy)         #s振幅反射率
    pRj=prj*(prj.conjugate())               #p能量反射率
    sRj=srj*(srj.conjugate())                #s能量反射率
    R=(pRj+sRj)/2
    return R.real





n,d,o,k=inputParameter()            # 输入层数k、折射率n、厚度d、入射角度o
Wavelength=wavelengthInput()     #要计算波长范围输入
reflectivity=[]           #能量反射率列表，空
#能量反射率计算
for wavelength in Wavelength:
    cos_ar_k=cos_RefractionAngle_k(o,n)      # cosφk的计算
    xwhd=PhaseThickness_k(d,cos_ar_k,wavelength,k,n)   #相位厚度
    PYK=Pyk(n,cos_ar_k,k)                     # p分量的光学导纳
    SYK=Syk(n,cos_ar_k,k)                    # s分量的光学导纳
    BCp=CharacteristicMatrix(PYK,xwhd,k)     # p分量矩阵[BC]计算
    BCs=CharacteristicMatrix(SYK,xwhd,k)     # s分量矩阵[BC]计算
    R=reflex(BCp,BCs,PYK[0],SYK[0])              #衬底和膜系组合的能量反射率
    reflectivity.append(R*100)
print(reflectivity)
plt.plot(Wavelength,reflectivity,'r')
plt.xlabel('Wavelength (nm)')      # 设置 x轴的名称
plt.ylabel('Reflectivity (%)')               # 设置 y轴的名称
plt.xlim(380,780)                # 设置x轴的坐标的范围
#plt.ylim(1,5)                  # 设置y轴的坐标的范围
plt.show()

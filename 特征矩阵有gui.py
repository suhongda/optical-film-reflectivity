# -*-coding:utf-8-*-
from tkinter import *
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import datetime


class main:

    def __init__(self):
        self.root = Tk()  # 实例化
        self.root.title('特征矩阵法')
        self.root.geometry('800x500')
        self.root.resizable(width=True , height=True)       # 可调整大小
        self.frm = Frame(self.root)      # 框架
        self.frm_top = Frame(self.frm)
        self.frm_mxcs = Frame(self.frm_top)

        Label(self.frm_mxcs, text='项目:', width=10, height=2,font=('微软雅黑',14),fg='red').grid(row=0, column=0)
        self.var_mxcs11 = StringVar()  # 接收输入
        self.entry_mxcs11 = Entry(self.frm_mxcs, textvariable=self.var_mxcs11, width=8,font=('微软雅黑',14),bd=0).grid(row=0, column=1)

        self.label_mxcs = Label(self.frm_mxcs,text='膜系层数:',width=10).grid(row=1, column=0)
        self.var_mxcs = StringVar()  # 接收输入
        self.entry_mxcs = Entry(self.frm_mxcs, textvariable=self.var_mxcs, width=13,bd=0).grid(row=1, column=1)        # 输入膜层
        self.button_mxcs = Button(self.frm_mxcs, text='确定', command=self.click,width=6,
                                                bd=0,bg='PaleGreen',activeforeground='yellow',activebackground='red').grid(row=1, column=2,padx=20)
        self.frm_mxcs.grid(row=1, column=0)
        #入射角
        self.frm_rsj = Frame(self.frm_top)
        self.label_rsj = Label(self.frm_rsj, text='入射角').grid(row=1, column=0, padx=5)
        self.var_rsj = DoubleVar(value=0)  # 接收入射角输入
        self.entry_rsj = Entry(self.frm_rsj, textvariable=self.var_rsj, width=5,bd=0).grid(row=1, column=1)
        self.frm_rsj.grid(row=1, column=1)
        # 波长范围
        self.frm_wave_range = Frame(self.frm_top)
        self.var_max = DoubleVar(value=780)  # 上限
        self.var_min = DoubleVar(value=380)  # 下限
        self.var_spacing = DoubleVar(value=1)  # 波长间隔
        self.label_wave_range = Label(self.frm_wave_range, text='波长范围:',width=15).grid(row=1, column=0)
        self.label_max = Label(self.frm_wave_range, text='上限:').grid(row=1, column=3)
        self.entry_range_max = Entry(self.frm_wave_range, textvariable=self.var_max, width=5, bd=0).grid(row=1, column=4)
        self.label_min = Label(self.frm_wave_range, text='下限:').grid(row=1, column=1)
        self.entry_range_min = Entry(self.frm_wave_range, textvariable=self.var_min, width=5, bd=0).grid(row=1, column=2)
        self.label_spacing = Label(self.frm_wave_range, text='波长间隔:').grid(row=1, column=5)
        self.entry_range_label_spacing = Entry(self.frm_wave_range, textvariable=self.var_spacing, width=5, bd=0).grid(row=1, column=6)
        self.frm_wave_range.grid(row=1, column=2, padx=20)
        self.frm_top.grid(row=1, sticky='w')

        # 底部
        self.frm_bottom = Frame(self.frm)
        # 特征矩阵
        self.frm_bottom_matrix = Frame(self.frm_bottom)
        # 画布
        self.canvas_left = Canvas(self.frm_bottom_matrix)
        # 动态生成的输入框,在画布里面
        self.frm_characteristic_matrix = Frame(self.canvas_left)
        Label(self.frm_characteristic_matrix, text='请先输入膜系层数，然后确认',font=('微软雅黑',19),fg='red').pack()
        # frm_characteristic_matrix绑定跟随画布滚动条滚动
        self.frm_characteristic_matrix.bind('<Configure>', self.scrollbar_event)
        # 滚动条
        self.scrollbar = Scrollbar(self.frm_bottom_matrix)
        self.scrollbar.pack(side=RIGHT, fill=Y)
        self.scrollbar['command'] = self.canvas_left.yview
        self.canvas_left.config(width=350)
        self.canvas_left.config(yscrollcommand=self.scrollbar.set)
        self.canvas_left.create_window((0, 0), window=self.frm_characteristic_matrix, anchor='nw')
        self.canvas_left.pack(side=LEFT, fill=BOTH)
        self.frm_bottom_matrix.pack(side=LEFT)
        # 折线图
        self.frm_bottom_chart = Frame(self.frm_bottom, height=400)
        # 新建画布
        # 在前面得到的子图上绘图
        # 将绘制的图形显示到tkinter:创建属于root的canvas画布,并将图f置于画布上
        self.figure, self.subplot = self.create_chart()
        self.canvas_chart = FigureCanvasTkAgg(self.figure, master=self.frm_bottom_chart)
        self.canvas_chart.draw()  # 注意show方法已经过时了,这里改用draw
        self.canvas_chart.get_tk_widget().pack(side=TOP,  # 上对齐
                                    fill=BOTH,  # 填充方式
                                    expand=YES)  # 随窗口大小调整而调整
        # 将绘制的图形显示到tkinter:创建画布,并将图f置于画布上
        self.frm_bottom_chart.pack(side=RIGHT)
        self.frm_bottom.grid(row=3, sticky='w')
        self.frm.pack()
        self.root.mainloop()

    def create_chart(self):
        f = Figure(figsize=(3, 3), dpi=145)
        f.subplots_adjust(left=0.2, bottom=0.2)
        a = f.add_subplot(111)  # 添加子图:1行1列第1个
        a.spines['bottom'].set_linewidth(0.5)  # 设置底部坐标轴的粗细
        a.spines['left'].set_linewidth(0.5) # 设置左边坐标轴的粗细
        a.spines['right'].set_linewidth(0.5)  # 设置右边坐标轴的粗细
        a.spines['top'].set_linewidth(0.5) # 设置上部坐标轴的粗细
        a.grid(linestyle='--',linewidth=0.4,which='major')        # 主刻度网格线
        a.grid(linestyle='--', linewidth=0.4, which='minor')      # 次刻度网格线
        a.tick_params(labelsize=7)      # 刻度的字体大小
        a.xaxis.set_minor_locator(MultipleLocator(20))    # 设置x轴次刻度间隔
        a.yaxis.set_minor_locator(MultipleLocator(5))   # 设置y轴次刻度间隔
        labels = a.get_xticklabels() + a.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]          # 刻度的字体
        a.set_xlabel('Wavelength (nm)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 x轴的名称
        a.set_ylabel('Reflectivity (%)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 y轴的名称
        a.set_xlim(380, 780)  # 设置x轴的坐标的范围
        a.set_ylim(0,100)                  # 设置y轴的坐标的范围
        return f,a

    def scrollbar_event(self, event):
        self.canvas_left.configure(scrollregion=self.canvas_left.bbox("all"))

    def click(self):
        self.mxcs = self.var_mxcs.get()
        self.gen()

    def gen(self):
        for widget in self.frm_characteristic_matrix.winfo_children():
            widget.destroy()
        mxcs_int = int(self.mxcs)
        self.n = [0 for x in range(mxcs_int+2)]
        self.d = [0 for x in range(mxcs_int+2)]
        for i in range(mxcs_int+2):
            if i == 0:
                label_level_1 = Label(self.frm_characteristic_matrix, text='入射介质n0：').grid(row=i, column=0,pady=2)
                self.n[i] = DoubleVar(value=1)
                entry_rsjz = Entry(self.frm_characteristic_matrix, textvariable=self.n[0], width=8,bd=0).grid(row=i, column=1)
            elif i == (mxcs_int+1):
                label_level_1 = Label(self.frm_characteristic_matrix, text='衬底n' + str(mxcs_int + 1) + '：').grid( row=mxcs_int + 1, column=0)
                self.n[i] = DoubleVar(value=1.52)
                entry_cd = Entry(self.frm_characteristic_matrix, textvariable=self.n[i], width=8,bd=0).grid(row=mxcs_int + 1, column=1)
                button_confirm = Button(self.frm_characteristic_matrix, text='确定', command=self.confirm,width=8,
                                                        bd=0,bg='PaleGreen',activeforeground='yellow',activebackground='red').grid(row=mxcs_int + 1, column=3)
            else:
                Label(self.frm_characteristic_matrix, text='第' + str(i) + '层n' + str(i) + ':').grid(row=i, column=0)
                self.n[i] = DoubleVar(value=1.45)
                entry_n = Entry(self.frm_characteristic_matrix, textvariable=self.n[i], width=8,bd=0).grid(row=i, column=1)
                Label(self.frm_characteristic_matrix,width=2).grid(row=i, column=2)
                Label(self.frm_characteristic_matrix, text='厚度d' + str(i) + ':').grid(row=i, column=3)
                self.d[i] = DoubleVar(value=100)
                entry_d = Entry(self.frm_characteristic_matrix, textvariable=self.d[i], width=7,bd=0).grid(row=i, column=4)
                Label(self.frm_characteristic_matrix, text='nm' ).grid(row=i, column=5)

    def draw_char(self, k, o, n, d, Wavelength):
        # 生成用于绘图的数据
        self.reflectivity = []  # 能量反射率列表，空
        # 能量反射率计算
        for wavelength in Wavelength:
            cos_ar_k = self.cos_RefractionAngle_k(o, n)  # cosφk的计算
            xwhd = self.PhaseThickness_k(d, cos_ar_k, wavelength, k, n)  # 相位厚度
            PYK = self.Pyk(n, cos_ar_k, k)  # p分量的光学导纳
            SYK = self.Syk(n, cos_ar_k, k)  # s分量的光学导纳
            BCp = self.CharacteristicMatrix(PYK, xwhd, k)  # p分量矩阵[BC]计算
            BCs = self.CharacteristicMatrix(SYK, xwhd, k)  # s分量矩阵[BC]计算
            R = self.reflex(BCp, BCs, PYK[0], SYK[0])  # 衬底和膜系组合的能量反射率
            self.reflectivity.append(R * 100)
        print('能量反射率：' + str(self.reflectivity))
        self.subplot.cla()
        self.subplot.plot(Wavelength, self.reflectivity, 'r',linewidth=0.5,alpha=1)
        self.subplot.spines['bottom'].set_linewidth(0.5)  # 设置底部坐标轴的粗细
        self.subplot.spines['left'].set_linewidth(0.5) # 设置左边坐标轴的粗细
        self.subplot.spines['right'].set_linewidth(0.5)  # 设置右边坐标轴的粗细
        self.subplot.spines['top'].set_linewidth(0.5) # 设置上部坐标轴的粗细
        self.subplot.xaxis.set_minor_locator(MultipleLocator(20))  # 设置x轴次刻度间隔
        self.subplot.yaxis.set_minor_locator(MultipleLocator(1))  # 设置y轴次刻度间隔
        self.subplot.grid(linestyle='--',linewidth=0.4,which='major')        # 网格线
        self.subplot.grid(linestyle='--', linewidth=0.4, which='minor')
        self.subplot.tick_params(labelsize=7)#
        labels = self.subplot.get_xticklabels() + self.subplot.get_yticklabels()#
        [label.set_fontname('Times New Roman') for label in labels]#
        self.subplot.set_xlabel('Wavelength (nm)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 x轴的名称
        self.subplot.set_ylabel('Reflectivity (%)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 y轴的名称
        # self.subplot.set_ylim(0:100)  # 设置y轴的坐标的范围
        self.canvas_chart.draw()




    def confirm(self):
        # 所有输入的数据都没有进行校验，当输入的数据错误是可能会出错
        k = int(self.var_mxcs.get())  # 膜系层数
        o = 2 * np.pi / 360 * float(self.var_rsj.get())  # 入射角度，入射角度（0°≤θ≤90°）
        n = [0] * (k+2)  # 初始化大小为k+2，全为0,多了入射介质，衬底
        n[0] = (float(self.n[0].get()))  # 先设置好 入射介质的折射率n0
        for i in range(1, k+2):
            n[i] = float(self.n[i].get())

        d1 = [0 for x in range(int(self.var_mxcs.get())+2)]   # 初始化为k+1 大小，为了使 n，d 的下标对应
        # d0 不知道存不存在，先设置为0吧
        d1[0] = 0
        for i in range(1, k + 1):
            d1[i] = float(self.d[i].get())
        # for i, j in self.n, self.d:
        #     if i.get() is None or j.get() is None:
        #         # 存在输入框为空，不能进行下去，不做过多处理
        #         return
        wave_length_min = self.var_min.get()
        wave_length_max = self.var_max.get()
        wave_length_spacing = self.var_spacing.get()
        # 范围
        Wavelength = self.wavelengthCaculate(wave_length_min, wave_length_max, wave_length_spacing)
        # 处理数据，没有对输入的数据进行验证，当输入的数据错误是可能会出错
        # 将处理好的数据当成参数，调用函数进行画图
        self.draw_char(k, o, n, d1, Wavelength) 
        # 输出图片
        fig = plt.figure()
        m = fig.add_subplot(111)
        plt.plot(Wavelength, self.reflectivity, 'r',linewidth=0.7,alpha=1)
        m.spines['bottom'].set_linewidth(0.5)  # 设置底部坐标轴的粗细
        m.spines['left'].set_linewidth(0.5) # 设置左边坐标轴的粗细
        m.spines['right'].set_linewidth(0.5)  # 设置右边坐标轴的粗细
        m.spines['top'].set_linewidth(0.5) # 设置上部坐标轴的粗细
        m.xaxis.set_major_locator(MultipleLocator(60))
        m.xaxis.set_minor_locator(MultipleLocator(20))  # 设置x轴次刻度间隔
        m.yaxis.set_minor_locator(MultipleLocator(1))  # 设置y轴次刻度间隔
        m.grid(linestyle='--',linewidth=0.4,which='major')        # 网格线
        m.grid(linestyle='--', linewidth=0.4, which='minor')
        m.tick_params(labelsize=7)#
        labels = m.get_xticklabels() + m.get_yticklabels()#
        [label.set_fontname('Times New Roman') for label in labels]#
        m.set_xlabel('Wavelength (nm)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 x轴的名称
        m.set_ylabel('Reflectivity (%)', fontdict={'family' : 'Times New Roman','size': 8})  # 设置 y轴的名
        time1 = datetime.datetime.now().strftime('%m.%d-%H%M%S.%f')
        plt.savefig(self.var_mxcs11.get()+'-'+time1+'.png', dpi=300)
        # 输出txt文件
        content1='项目:'+self.var_mxcs11.get()
        content2='时间:'+datetime.datetime.now().strftime('%Y.%m.%d-%H:%M:%S.%f')
        content3='层数:'+str(k)+'        '+'入射角:'+str(o)
        content4=''
        content5 = ''
        for i in range(int(self.mxcs)+2):
            if i==0  :
                co1=str(i)+' '+str(n[i])+'   --'+'\n'
                content4=content4+co1
            elif i==int(self.mxcs)+1:
                co1 = str(i) + ' ' + str(n[i]) + '   --' + '\n'
                content4 = content4 + co1
            else:
                co1 = str(i) + ' ' + str(n[i]) +' '+str(d1[i])+ '\n'
                content4 = content4 + co1
        for i in range(int( (wave_length_max-wave_length_min)/wave_length_spacing)):
            co2=str(Wavelength[i]) + ' ' + str(self.reflectivity[i])+ '\n'
            content5 = content5 + co2
        # 写入文件
        with open(self.var_mxcs11.get()+'-'+time1 + '.txt', 'w') as f:
            f.write(content1+'\n'
                    +content2+'\n'
                    +content3+'\n\n'
                    +'层  折射率  厚度\n'
                    +content4+'\n'
                    +'波长  反射率\n'
                    +content5)



    # 要计算波长范围输入
    def wavelengthCaculate(self, wmin, wmax, s):  # 要计算波长范围输入
        wavelengthMin = int(wmin)
        spacing = float((s))
        wavelengthMax = int(wmax)
        Wavelength = np.arange(wavelengthMin, wavelengthMax, spacing)  # 产生波长序列
        return Wavelength

    # 以下复制过来的
    # cosφk的计算
    def cos_RefractionAngle_k(self,o, n):  # cosφk的计算
        cos_RA = []
        sin_RAe_0 = np.sin(o)  # sinφ0 的计算
        for i in n:
            sin_RAe_k = n[0] * sin_RAe_0 / i
            cos_RefractionAngle = (1 - sin_RAe_k ** 2) ** 0.5  # cosφk=(1-(n0*sinφ0 / nk)²)^0.5
            cos_RA.append(cos_RefractionAngle)
        return cos_RA

    # 相位厚度
    def PhaseThickness_k(self,d, cos_ar_k, wavelength, k, n):  # 相位厚度
        ptk = []
        kk = 0
        while kk <= k + 1:
            if kk == 0:
                ptk.append(0)  # 入射介质相位厚度
            else:
                PT_k = 2 * np.pi * cos_ar_k[kk] * n[kk] * d[kk] / wavelength  # 2π*dk*cosφk / λ
                ptk.append(PT_k)
            kk = kk + 1
        return ptk

    # p分量的光学导纳
    def Pyk(self,n, cos_ar_k, k):  # p分量的光学导纳
        hdsj = 0
        pyk = []
        while hdsj <= k + 1:
            pyk.append(n[hdsj] / cos_ar_k[hdsj])  # ηk=nk / cosφk
            hdsj = hdsj + 1
        return pyk

    # s分量的光学导纳
    def Syk(self,n, cos_ar_k, k):  # s分量的光学导纳
        hdsj1 = 0
        syk = []
        while hdsj1 <= k + 1:
            syk.append(n[hdsj1] * cos_ar_k[hdsj1])  # ηk=nk*cosφk
            hdsj1 = hdsj1 + 1
        return syk

    # 矩阵[BC]计算
    def CharacteristicMatrix(self,YK, xwhd, k):  # 矩阵[BC]
        bc_2 = [1, YK[-1]]
        gdsh = 0
        bc_1 = []
        while gdsh <= k:  # 各层特征矩阵连乘
            if gdsh == 0:
                a0 = np.identity(2)  # 创建单位矩阵
                bc_1.append(a0)
            else:
                m11 = np.cos(xwhd[gdsh])
                m22 = np.cos(xwhd[gdsh])  # m11,m12,m21,m22特征矩阵的各元素
                m12 = 1j * np.sin(xwhd[gdsh]) / YK[gdsh]
                m21 = 1j * np.sin(xwhd[gdsh]) * YK[gdsh]
                a1 = np.array([[m11, m12], [m21, m22]])  # 创建特征矩阵
                a2 = np.dot(a1, bc_1[-1])  # 前一个特征矩阵与后一个特征矩阵
                bc_1.append(a2)
            gdsh = gdsh + 1
        bc_11 = bc_1[-1]
        bc = np.dot(bc_11, bc_2)  # 得到BC矩阵
        return bc

    # 衬底和膜系组合的能量反射率
    def reflex(self,BCp, BCs, PYK_0, SYK_0):  # 衬底和膜系组合的能量反射率
        py = BCp[1] / BCp[0]  # p的等效导纳
        sy = BCs[1] / BCs[0]  # s的等效导纳
        prj = (PYK_0 - py) / (PYK_0 + py)  # p振幅反射率
        srj = (SYK_0 - sy) / (SYK_0 + sy)  # s振幅反射率
        pRj = prj * (prj.conjugate())  # p能量反射率
        sRj = srj * (srj.conjugate())  # s能量反射率
        R = (pRj + sRj) / 2
        return R.real


if __name__ == '__main__':
    main()

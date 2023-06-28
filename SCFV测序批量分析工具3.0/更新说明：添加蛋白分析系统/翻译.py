from Bio.Seq import Seq
from Bio import SeqIO
import re
import os
import tkinter
from tkinter import Tk, StringVar, filedialog, IntVar, Text
import vl_align
from tkinter.ttk import *


class Application_ui(Frame):
    #这个类仅实现界面生成功能，具体事件处理代码在子类Application中。

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master.title('SCFV分析软件3.0')
        self.master.geometry('772x516')
        # self.master.iconphoto(False, tkinter.PhotoImage(file='bbb.png'))

        self.createWidgets()

    def createWidgets(self):
        self.style = Style()
        self.value1 = IntVar()
        self.value2 = IntVar()
        self.value3 = IntVar(value=None)

        tabControl = Notebook(self.master)
        self.tab1 = Frame(tabControl)
        tabControl.add(self.tab1, text='     Scfv 序列分析     ')

        self.tab2 = Frame(tabControl)
        tabControl.add(self.tab2, text='      抗体蛋白分析      ')

        self.style.configure('Frame4.TLabelframe', font=('宋体', 9))
        self.Frame4 = LabelFrame(self.tab1, text='结果输出', style='Frame4.TLabelframe')
        self.Frame4.place(relx=0.518, rely=0.047, relwidth=0.416, relheight=0.886)

        # self.style.configure('Command1.TButton', font=('宋体', 9))
        # self.Command1 = Button(self.top, text='输出', command=self.Command2_Cmd, style='Command1.TButton')
        # self.Command1.place(relx=0.249, rely=0.868, relwidth=0.136, relheight=0.064)

        self.style.configure('Command1.TButton', font=('宋体', 9))
        self.Command2 = Button(self.tab1, text='分析全长', command=self.Command1_Cmd, style='Command1.TButton')
        self.Command2.place(relx=0.035, rely=0.868, relwidth=0.136, relheight=0.064)

        self.style.configure('Command1.TButton', font=('宋体', 9))
        self.Command2 = Button(self.tab1, text='分别分析', command=self.Command2_Cmd, style='Command1.TButton')
        self.Command2.place(relx=0.175, rely=0.868, relwidth=0.136, relheight=0.064)

        self.style.configure('Command1.TButton', font=('宋体', 9))
        self.Command2 = Button(self.tab1, text='错误信息', command=self.Command3_Cmd, style='Command1.TButton')
        self.Command2.place(relx=0.315, rely=0.868, relwidth=0.136, relheight=0.064)

        self.style.configure('Frame3.TLabelframe', font=('宋体',9))
        self.Frame3 = LabelFrame(self.tab1, text='请选择命名方式', style='Frame3.TLabelframe')
        self.Frame3.place(relx=0.041, rely=0.589, relwidth=0.405, relheight=0.25)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option3 = Radiobutton(self.Frame3, text='输出全部名称', variable=self.value1, value=1, style='Option1.TRadiobutton')
        self.Option3.place(relx=0.128, rely=0.301, relwidth=0.588, relheight=0.548)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option4 = Radiobutton(self.Frame3, text='输出‘_’前名称', variable=self.value1, value=0, style='Option1.TRadiobutton')
        self.Option4.place(relx=0.128, rely=0.153, relwidth=0.588, relheight=0.348)

        self.style.configure('Frame2.TLabelframe', font=('宋体', 9))
        self.Frame2 = LabelFrame(self.tab1, text='请选择是否为拼接序列', style='Frame2.TLabelframe')
        self.Frame2.place(relx=0.041, rely=0.31, relwidth=0.405, relheight=0.25)

        self.style.configure('Frame1.TLabelframe', font=('宋体', 9))
        self.Frame1 = LabelFrame(self.tab1, text='请选择文件夹', style='Frame1.TLabelframe')
        self.Frame1.place(relx=0.041, rely=0.047, relwidth=0.405, relheight=0.234)

        self.s1 = tkinter.Scrollbar(self.Frame4, orient=tkinter.VERTICAL)
        self.s1.pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.s1.config()

        self.Text2Var = StringVar(value='Text2')
        self.Text2 = Text(self.Frame4, font=('Times New Roman', 10), yscrollcommand=self.s1.set, wrap=tkinter.NONE)
        self.Text2.place(relx=0.025, rely=0.035, relwidth=0.95, relheight=0.947)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option1 = Radiobutton(self.Frame2, text='不是拼接序列', variable=self.value2, value=1, style='Option1.TRadiobutton')
        self.Option1.place(relx=0.128, rely=0.496, relwidth=0.31, relheight=0.194)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option2 = Radiobutton(self.Frame2, text='是拼接序列', variable=self.value2, value=0, style='Option1.TRadiobutton')
        self.Option2.place(relx=0.128, rely=0.248, relwidth=0.31, relheight=0.194)

        self.select_path1 = StringVar()
        self.Text1Var = StringVar(value='Text1')
        self.Text1 = Entry(self.Frame1, textvariable=self.select_path1, font=('宋体', 9))
        self.Text1.place(relx=0.051, rely=0.363, relwidth=0.642, relheight=0.207)

        self.style.configure("filedialog", font=('宋体', 9))
        self.Command2 = Button(self.tab1, text='选择', command=self.select_file1)
        self.Command2.place(relx=0.341, rely=0.155, relwidth=0.078, relheight=0.045)


        self.style.configure('Frame1.TLabelframe', font=('宋体', 9))
        self.Frame6 = LabelFrame(self.tab2, text='请选择文件夹', style='Frame1.TLabelframe')
        self.Frame6.place(relx=0.041, rely=0.047, relwidth=0.405, relheight=0.234)

        self.select_path2 = StringVar()
        self.Text3Var = StringVar(value='Text3')
        self.Text3 = Entry(self.Frame6, textvariable=self.select_path2, font=('宋体', 9))
        self.Text3.place(relx=0.051, rely=0.363, relwidth=0.642, relheight=0.207)
        self.style.configure("filedialog", font=('宋体', 9))
        self.Command5 = Button(self.tab2, text='选择', command=self.select_file2)
        self.Command5.place(relx=0.341, rely=0.155, relwidth=0.078, relheight=0.045)


        self.style.configure('Frame2.TLabelframe', font=('宋体', 9))
        self.Frame5 = LabelFrame(self.tab2, text='请选择轻链或重链', style='Frame2.TLabelframe')
        self.Frame5.place(relx=0.041, rely=0.31, relwidth=0.405, relheight=0.25)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option5 = Radiobutton(self.Frame5, text='轻链抗体VL', variable=self.value3, value=2,
                                   style='Option1.TRadiobutton')
        self.Option5.place(relx=0.128, rely=0.496, relwidth=0.31, relheight=0.194)

        self.style.configure('Option1.TRadiobutton', font=('宋体', 9))
        self.Option6 = Radiobutton(self.Frame5, text='重链抗体VH', variable=self.value3, value=1,
                                   style='Option1.TRadiobutton')
        self.Option6.place(relx=0.128, rely=0.248, relwidth=0.31, relheight=0.194)


        self.style.configure('Frame4.TLabelframe', font=('宋体', 9))
        self.Frame7 = LabelFrame(self.tab2, text='结果输出', style='Frame4.TLabelframe')
        self.Frame7.place(relx=0.518, rely=0.047, relwidth=0.416, relheight=0.886)

        self.Text4Var = StringVar(value='Text2')
        self.Text4 = Text(self.Frame7, font=('Times New Roman', 10), yscrollcommand=self.s1.set, wrap=tkinter.NONE)
        self.Text4.place(relx=0.025, rely=0.035, relwidth=0.95, relheight=0.947)

        self.style.configure('Command1.TButton', font=('宋体', 9))
        self.Command2 = Button(self.tab2, text='分析', command=self.align, style='Command1.TButton')
        self.Command2.place(relx=0.175, rely=0.868, relwidth=0.136, relheight=0.064)

        tabControl.pack(expand=1, fill="both")

    def select_file1(self):
        selected_folder = filedialog.askdirectory()  # 使用askdirectory函数选择文件夹
        self.select_path1.set(selected_folder)
    def select_file2(self):
        selected_folder = filedialog.askdirectory()  # 使用askdirectory函数选择文件夹
        self.select_path2.set(selected_folder)


class Application(Application_ui):
    #这个类实现具体的事件处理回调函数。界面生成代码在Application_ui中。
    def __init__(self, master=None):
        Application_ui.__init__(self, master)
        self.message = []

    def Command1_Cmd(self):

        self.Text2.delete('1.0', tkinter.END)

        path = str(self.select_path1.get())
        pinjie = int(self.value2.get())
        out_name = int(self.value1.get())
        vl_count = 0
        vh_count = 0


        if not path:
            self.Text2.insert(tkinter.END, '请选择文件夹')
            return

        for i in InPut(path).listdir()[0]:
            if int(pinjie) == 0:
                # 拼接序列的文件名输出
                name = OutPut(path='', name=i).name_moth_1()

                if int(out_name) == 1:
                    me = Fenxi(i, path, pinjie, outname=i[:-4]).do_complay()
                    if me[0]:
                        self.Text2.insert(tkinter.END, me[0] + '\n')
                    if me[1]:
                        self.Text2.insert(tkinter.END, me[1])
                    if me[2] == 1:
                        vh_count += 1
                    if me[3] == 1:
                        vl_count += 1
                    if me[4]:
                        self.message.append(me[4])
                    if me[5]:
                        self.message.append(me[5])

                elif int(out_name) == 0:
                    try:
                        me = Fenxi(i, path, pinjie, outname=name[0]).do_complay()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0] + '\n')
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])
                    except:
                        self.Text2.insert(tkinter.END, '设置存在问题，按原名输出'+'\n')
                        me = Fenxi(i, path, pinjie, outname=name[0]).do_complay()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0] + '\n')
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])

            elif int(pinjie) == 1:
                # 非拼接序列的文件名输出
                name = OutPut(path='', name=i).name_moth_2()
                if int(out_name) == 1:
                    me = Fenxi(i, path, pinjie, outname=i[:-4]).do_complay()
                    if me[0]:
                        self.Text2.insert(tkinter.END, me[0] + '\n')
                    if me[1]:
                        self.Text2.insert(tkinter.END, me[1])
                    if me[2] == 1:
                        vh_count += 1
                    if me[3] == 1:
                        vl_count += 1
                    if me[4]:
                        self.message.append(me[4])
                    if me[5]:
                        self.message.append(me[5])
                elif int(out_name) == 0:
                    try:
                        me = Fenxi(i, path, pinjie, outname=name[0][0] + '_' + name[0][1]).do_complay()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0] + '\n')
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])
                    except:

                        self.Text2.insert(tkinter.END, '设置存在问题，按原名输出'+'\n')
                        me = Fenxi(i, path, pinjie, outname=i[:-4]).do_complay()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0] + '\n')
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])

        self.Text2.insert(tkinter.END, '\n' + '成功解析序列:'+str(vl_count)+'/'+str(len(InPut(path).listdir()[0])) + '\n')

    def Command2_Cmd(self, event=None):
        self.Text2.delete('1.0', tkinter.END)
        path = str(self.select_path1.get())
        pinjie = int(self.value2.get())
        out_name = int(self.value1.get())
        vl_count = 0
        vh_count = 0

        if not path:
            self.Text2.insert(tkinter.END, '请选择文件夹')
            return

        for i in InPut(path).listdir()[0]:
            if int(pinjie) == 0:
                # 拼接序列的文件名输出
                name = OutPut(path='', name=i).name_moth_1()

                if int(out_name) == 1:
                    me = Fenxi(i, path, pinjie, outname=i[:-4]).do_vhvl()
                    if me[0]:
                        self.Text2.insert(tkinter.END, me[0])
                    if me[1]:
                        self.Text2.insert(tkinter.END, me[1])
                    if me[2] == 1:
                        vh_count += 1
                    if me[3] == 1:
                        vl_count += 1
                    if me[4]:
                        self.message.append(me[4])
                    if me[5]:
                        self.message.append(me[5])

                elif int(out_name) == 0:
                    try:
                        me = Fenxi(i, path, pinjie, outname=name[0]).do_vhvl()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0])
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])
                    except:
                        self.Text2.insert(tkinter.END, '\n' + '设置存在问题，按原名输出')

                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0])
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                            return
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])

            elif int(pinjie) == 1:
                # 非拼接序列的文件名输出
                name = OutPut(path='', name=i).name_moth_2()
                if int(out_name) == 1:
                    me = Fenxi(i, path, pinjie, outname=i[:-4]).do_vhvl()
                    if me[0]:
                        self.Text2.insert(tkinter.END, me[0])
                    if me[1]:
                        self.Text2.insert(tkinter.END, me[1])
                    if me[2] == 1:
                        vh_count += 1
                    if me[3] == 1:
                        vl_count += 1
                    if me[4]:
                        self.message.append(me[4])
                    if me[5]:
                        self.message.append(me[5])
                elif int(out_name) == 0:
                    try:
                        me = Fenxi(i, path, pinjie, outname=name[0][0] + '_' + name[0][1]).do_vhvl()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0])
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])
                    except:

                        self.Text2.insert(tkinter.END, '\n' + '设置存在问题，按原名输出')
                        me = Fenxi(i, path, pinjie, outname=i[:-4]).do_vhvl()
                        if me[0]:
                            self.Text2.insert(tkinter.END, me[0])
                        if me[1]:
                            self.Text2.insert(tkinter.END, me[1])
                        if me[2] == 1:
                            vh_count += 1
                        if me[3] == 1:
                            vl_count += 1
                        if me[4]:
                            self.message.append(me[4])
                        if me[5]:
                            self.message.append(me[5])

        self.Text2.insert(tkinter.END, '\n' + '成功解析vh序列:' + str(vh_count) + '/' + str(
            len(InPut(path).listdir()[0])) + '\n' + '成功解析vl序列:' + str(vl_count) + '/' + str(
            len(InPut(path).listdir()[0])) + '\n')

    def Command3_Cmd(self, event=None):
        self.Text2.delete('1.0', tkinter.END)
        if not self.message:
            self.Text2.insert(tkinter.END, '恭喜你未查询到错误信息！'+'\n')
        else:
            self.Text2.insert(tkinter.END, '\n\n')
            for i in self.message:
                self.Text2.insert(tkinter.END, i + '\n')

    def align(self):
        path = str(self.select_path2.get())
        if not path:
            self.Text4.delete('1.0', tkinter.END)
            self.Text4.insert(tkinter.END, '请选择文件夹'+'\n')
            return

        name, path1 = vl_align.InPut(path).listdir()
        h_l = int(self.value3.get())
        print(h_l)
        if not h_l:
            self.Text4.delete('1.0', tkinter.END)
            self.Text4.insert(tkinter.END, '请选择轻链或重链' + '\n')
            return
        A = vl_align.Vl_Align(path, h_l)
        for i in name:
            A.read_pro(i)
        A.align()
        A.do()
        A.save_file(path)

        self.Text4.delete('1.0', tkinter.END)
        for i in A.message:
            self.Text4.insert(tkinter.END, i + '\n')


class InPut(object):  # 读取序列文件

    def __init__(self, __path):
        self.file_name = []
        self.__path = __path

    def listdir(self):
        # 读取.seq序列文件
        for __file in os.listdir(self.__path):
            if '.seq' in __file:
                __file_path = os.path.join(self.__path, __file)
                self.file_name.append(__file)

        # 读取.fasta文件
        # elif '.fasta' in __file:
        #     __file_path = os.path.join(self.__path, __file)
        #     print(__file_path)

        # print(self.file_name)
        return self.file_name, self.__path


class Fenxi(object):  # 分析序列

    def __init__(self, file_name, file_path, pinjie, outname):

        self.vh_messages = ''
        self.vl_messages = ''
        self.messages = ''
        self.pro = []
        self.my_seq = ''
        self.seq_set = {}
        self.out_name = outname
        self.scfv_seq = []
        self.vh_hid_pro = ReadTest().read_vh_pro_hid()
        self.vh_foot_pro = ReadTest().read_vh_pro_foot()
        self.vl_hid_pro = ReadTest().read_vl_pro_hid()
        self.vl_foot_pro = ReadTest().read_vl_pro_foot()
        self.vh_hid_seq = ''
        self.vh_foot_seq = ''
        self.vl_hid_seq = ''
        self.vl_foot_seq = ''
        self.complay_seq = ''
        self.file_name = file_name
        self.file_path = file_path
        self.pinjie = int(pinjie)
        self.error = ''
        self.vh_count = 0
        self.vl_count = 0

    def duqu(self):
        # 读取序列
        if self.pinjie == 0:  # 读取拼接序列
            self.my_seq = Seq(str(self.pinjieseq()))
        else:  # 读取非拼接序列
            self.my_seq = Seq(str(ReadTest.read_seq(self.file_name, self.file_path)))

        # 若序列反向则反向互补获得正向序列
        if 'GGCCGCCTG' in self.my_seq:
            self.my_seq = self.my_seq.reverse_complement()
        elif 'GGCCGGCCT' in self.my_seq:
            self.my_seq = self.my_seq.reverse_complement()

    def pinjieseq(self):
        # 读取拼接序列
        try:
            seq = SeqIO.read(self.file_path + '\\' + self.file_name, 'fasta').seq
            return seq
        except:
            seq = ReadTest.read_seq(self.file_name, self.file_path)
            self.error = '不是拼接序列，按非拼接序列输出:' + self.out_name + '\n'
            return seq

    def fanyi(self, ):
        # 序列翻译为蛋白序列
        for i in [0, 1, 2]:
            s = self.my_seq[i:]
            while len(s) % 3 != 0:
                s = s + 'A'
            else:
                try:
                    a = s.translate(stop_symbol='.')
                    self.pro.append(a)
                except:

                    self.error = '请检查序列是否对应选项'
                    return

    def vh_hid(self):  # 读取配制文件，匹配起始位点
        for j in self.pro:
            for i in self.vh_hid_pro:
                seq = re.findall(r'%s.*' % i, str(j), re.M | re.I | re.M)
                if seq:
                    self.vh_hid_seq = (str(seq[0]))

    def vh_foot(self):  # 读取配制文件，匹配结束位点
        if self.vh_hid_seq:
            for j in self.vh_foot_pro:
                seq = re.findall(r'.*%s' % j, self.vh_hid_seq)
                if seq:
                    self.vh_foot_seq = str(seq[0])
                    return
            else:
                self.vh_messages = 'vh序列错误请核对：' + self.out_name
                return 'vh序列错误请核对：' + self.out_name
        else:
            # print('未找到vh头部序列，请手动核对：'+self.out_name)
            self.vh_messages = '未找到vh序列请核对：' + self.out_name
            return '未找到vh序列请核对：' + self.out_name

    def vl_hid(self):
        for j in self.pro:
            for i in self.vl_hid_pro:
                seq = re.findall(r'%s.*' % i, str(j), re.M | re.I | re.M)

                if seq:
                    self.vl_hid_seq = (str(seq[0]))
                    return

    def vl_foot(self):  # 读取配制文件，匹配结束位点
        seq = ''
        if self.vl_hid_seq:
            for j in self.vl_foot_pro:
                seq = re.findall(r'.*%s' % j, self.vl_hid_seq)
                if seq:
                    self.vl_foot_seq = str(seq[0])
                    return
            else:
                self.vl_messages = 'vl 序列错误请核对：' + self.out_name
                # print('未找到vl尾部序列，请手动核对：'+self.out_name)
                return 'vl 序列错误请核对：' + self.out_name
        else:
            # print('未找到vl起始序列，请手动核对：'+self.out_name)
            self.vl_messages = '未找到vl序列请核对：' + self.out_name
            return

    def vl_foot_complay(self):  # 读取配制文件，匹配结束位点

        if not self.vh_hid_seq:
            # print('未找到vl起始序列，请手动核对：'+self.out_name)
            self.vl_messages = '未找到序列请核对：' + self.out_name
            return

        for j in self.vl_foot_pro:
            seq = re.findall(fr'.*{j}', str(self.vh_hid_seq), re.M | re.I)
            if seq:
                self.vl_foot_seq = str(seq[0])
                return

        if not self.vl_foot_seq:
            self.vl_messages = '未找到尾部序列请核对：' + self.out_name
            # print('未找到vl尾部序列，请手动核对：'+self.out_name)
            return '未找到尾部序列请核对：' + self.out_name

    def shaixuan(self):
        # 筛选出正确序列
        if self.vl_foot_seq:
            if 'GGGGSGGGGSGGGGS' in self.vl_foot_seq:
                self.complay_seq = self.vl_foot_seq
                return 1
            elif 'GGSGGSGGSGGSGGSGGSGGS' in self.vl_foot_seq:  # 筛选阳参序列
                self.complay_seq = self.vl_foot_seq
            else:
                self.messages = 'linker序列错误：' + self.out_name
        elif self.vh_hid_seq:
            if not 'GGGGSGGGGSGGGGS' in self.vh_hid_seq:
                self.messages = '序列错误：' + self.out_name
            else:
                self.messages = '尾部序列错误：' + self.out_name
        else:
            self.messages = '未知序列：' + self.out_name

    def qishihid(self):
        # 输出起始代码(测试用）
        vh_hid_lis = []
        vh_hid_count = []

        for seq in self.scfv_seq:
            my_seq_pro = seq[0:5]
            vh_hid_lis.append(str(my_seq_pro))
            # print(vh_hid_lis)

        for i in vh_hid_lis:
            seq_count = vh_hid_lis.count(i)
            if [i, seq_count] not in vh_hid_count:
                self.seq_set[i] = seq_count
        print('vh_hid:', self.seq_set)

    def vh_write_pro(self, vh_pro):
        if vh_pro:
            try:
                os.makedirs('%s\\%s scFv蛋白序列\\vh 蛋白序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv蛋白序列\\vh 蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vh_pro)
                self.vh_count = 1
            except:
                with open('%s\\%s scFv蛋白序列\\vh 蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vh_pro)
                self.vh_count = 1

    def vl_write_pro(self, vl_pro):
        if vl_pro:
            try:
                os.makedirs('%s\\%s scFv蛋白序列\\vl 蛋白序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv蛋白序列\\vl 蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vl_pro)
                self.vl_count = 1
            except:
                with open('%s\\%s scFv蛋白序列\\vl 蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vl_pro)
                self.vl_count = 1

    def vh_write_seq(self, vh_seq):
        if vh_seq:
            try:
                os.makedirs('%s\\%s scFv基因序列\\vh 基因序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv基因序列\\vh 基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vh_seq)
            except:
                with open('%s\\%s scFv基因序列\\vh 基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vh_seq)

    def vl_write_seq(self, vl_seq):
        if vl_seq:
            try:
                os.makedirs('%s\\%s scFv基因序列\\vl 基因序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv基因序列\\vl 基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vl_seq)
            except:
                with open('%s\\%s scFv基因序列\\vl 基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(vl_seq)

    def complay_write_pro(self, complay_pro):
        if complay_pro:
            try:
                os.makedirs('%s\\%s scFv蛋白序列\\蛋白序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv蛋白序列\\蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(complay_pro)
                self.vl_count = 1
            except:
                with open('%s\\%s scFv蛋白序列\\蛋白序列\\%s.pro' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(complay_pro)
                self.vl_count = 1

    def complay_write_seq(self, complay_seq):
        if complay_seq:
            try:
                os.makedirs('%s\\%s scFv基因序列\\基因序列' % (self.file_path, self.file_name[:5]))
                with open('%s\\%s scFv基因序列\\基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(complay_seq)
            except:
                with open('%s\\%s scFv基因序列\\基因序列\\%s.seq' % (self.file_path, self.file_name[:5], self.out_name), 'w+') as f:
                    f.write(complay_seq)

    def vh_jiequseq(self):
        # 截取序列，获得正确序列得基因序列
        vh_good_seq = re.findall(r'CAGGCGGCC(.*)GCCCCATCG', str(self.my_seq), )
        if vh_good_seq:
            return vh_good_seq[0]+'GCCCCATCG'
        else:
            vh_good_seq = re.findall(r'CAGGCGGCC(.*)TGGGGCCTT', str(self.my_seq), )
            if vh_good_seq:
                return vh_good_seq[0]+'TGGGGCCTT'
            else:
                self.vh_messages = 'vh截取序列错误，请手动核对：' + self.out_name
                return

                # return 'vh序列错误请核对：' + self.out_name

    def vl_jiequseq(self):
        # 截取序列，获得正确序列得基因序列
        vl_good_seq = re.findall(r'GGCGGATCG(.*)GGCCAGG', str(self.my_seq), )
        if vl_good_seq:
            return vl_good_seq[0]
        else:
            self.vl_messages = 'vl截取序列错误，请手动核对：'+self.out_name
            return

    def complay_jiequseq(self):
        vl_good_seq = re.findall(r'CAGGCGGCC(.*)GGCCAGG', str(self.my_seq), )
        if vl_good_seq:
            return vl_good_seq[0]
        else:
            self.vl_messages = '截取序列错误，请手动核对：' + self.out_name
            return

    def do_vhvl(self):
        # 执行程序
        self.duqu()
        self.fanyi()
        self.vh_hid()
        self.vh_foot()
        self.vl_hid()
        self.vl_foot()
        if self.vh_foot_seq:
            self.vh_write_pro(self.vh_foot_seq)
            self.vh_write_seq(self.vh_jiequseq())
        if self.vl_foot_seq:
            self.vl_write_pro(self.vl_foot_seq)
            self.vl_write_seq(self.vl_jiequseq())
        return self.messages, self.error, self.vh_count, self.vl_count, self.vh_messages, self.vl_messages

    def do_complay(self):
        self.duqu()
        self.fanyi()
        self.vh_hid()
        self.vl_foot_complay()
        self.shaixuan()
        self.complay_write_pro(self.complay_seq)
        self.complay_write_seq(self.complay_jiequseq())
        return self.messages, self.error, self.vh_count, self.vl_count, self.vh_messages, self.vl_messages


class ReadTest(object):

    def read_seq(filename, filepath):

        with open(filepath + '\\' + filename, 'r+') as f:
            r_seq = f.read().replace('\t', '').replace(' ', '')
            if '^' in r_seq:
                true_seq =r_seq.replace(r_seq[0:r_seq.find('^\n')], '').replace('^\n', '')
            elif r_seq[-1] == '/' and 'ORIGIN' in r_seq:
                pretrue_seq = str(r_seq[r_seq.find('ORIGIN') + 6:-1].strip('/'))
                true_seq = re.sub(r'\d', '', pretrue_seq)
            else:
                true_seq = r_seq

            return true_seq

    def read_vh_pro_hid(self):
        pro_hid_list = []
        try:
            with open('配制文件.txt', 'r', encoding='utf-8') as f:
                a = re.findall(r'vh头部序列:\n(.*)', str(f.read()), re.M | re.I)
                pro_hid = re.findall(r'([a-z]+),', str(a), re.M | re.I)
                for i in pro_hid:
                    pro_hid_list.append(i)
                # print(pro_hid_list)
            return pro_hid_list
        except:
            print('未查询到配制文件，请将配制文件放置同一目录下。')
            return ['未查询到配制文件，请将配制文件放置同一目录下。']


    def read_vh_pro_foot(self):
        pro_foot_list = []
        with open('配制文件.txt', 'r', encoding='utf-8') as f:
            a = re.findall(r'vh尾部序列:\n(.*)', str(f.read()), re.M | re.I)
            pro_foot = re.findall(r'([a-z]+),', str(a), re.M | re.I)
            for i in pro_foot:
                pro_foot_list.append(i)
            # print(pro_foot_list)
        return pro_foot_list

    def read_vl_pro_hid(self):
        pro_hid_list = []
        with open('配制文件.txt', 'r', encoding='utf-8') as f:
            a = re.findall(r'vl头部序列:\n(.*)', str(f.read()), re.M | re.I)
            pro_hid = re.findall(r'([a-z]+),', str(a), re.M | re.I)
            for i in pro_hid:
                pro_hid_list.append(i)
            # print(pro_hid_list)
        return pro_hid_list



    def read_vl_pro_foot(self):
        pro_foot_list = []
        with open('配制文件.txt', 'r', encoding='utf-8') as f:
            a = re.findall(r'vl尾部序列:\n(.*)', str(f.read()), re.M | re.I)
            pro_foot = re.findall(r'([a-z]+),', str(a), re.M | re.I)
            for i in pro_foot:
                pro_foot_list.append(i)
            # print(pro_foot_list)
        return pro_foot_list


class OutPut(object):
    def __init__(self, path, name):
        self.name_moth = 0
        self.name = name
        self.path_moth = 0
        self.path = path

    def name_moth_1(self):  # 获得_前名称
        name = re.findall(r'(.*?)_.*?', str(self.name), re.M | re.I)
        return name

    def name_moth_2(self):
        name = re.findall(r'(.*?)_(.*?)_', str(self.name), re.M | re.I)
        return name


if __name__ == '__main__':

    top = Tk()
    Application(top).mainloop()
    try:
        top.destroy()
    except: pass


#coding:utf-8
"""
Created on Tue Nov 20 12:41:40 2018
@author: yujin.wang
"""
from Tkinter import *
import tkMessageBox
import tushare as ts

root = Tk()
root.title("hello world")
root.geometry('300x300')#是x 不是*


l1 = Label(root, text="xls名：")
l1.pack()#这里的side可以赋值为LEFT  RTGHT TOP  BOTTOM
xls_text = StringVar()
xls = Entry(root, textvariable = xls_text)
xls_text.set(" ")
xls.pack()


l2 = Label(root, text="sheet名：")
l2.pack()#这里的side可以赋值为LEFT  RTGHT TOP  BOTTOM
sheet_text = StringVar()
sheet = Entry(root, textvariable = sheet_text)
sheet_text.set(" ")
sheet.pack()


l3 = Label(root, text="循环次数：")
l3.pack()#这里的side可以赋值为LEFT  RTGHT TOP  BOTTOM
loop_text = StringVar()
loop = Entry(root, textvariable = loop_text)
loop_text.set(" ")
loop.pack()




l4 = Label(root, text="休眠时间：")
l4.pack()#这里的side可以赋值为LEFT  RTGHT TOP  BOTTOM
sleep_text = StringVar()
sleep = Entry(root, textvariable = sleep_text)
sleep_text.set(" ")
sleep.pack()


def on_click():
    x = xls_text.get()
    s = sheet_text.get()
    l = loop_text.get()
    sl = sleep_text.get()
    string = str("xls名：%s sheet名：%s 循环次数：%s 休眠时间：%s " %(x, s, l, sl))
    # print string
    df = str(ts.get_tick_data('600848',date='2018-12-12',src='tt'))
    tkMessageBox.showinfo(title='aaa', message = df)
Button(root, text="press", command = on_click).pack()
root.mainloop()

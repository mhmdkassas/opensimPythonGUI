from tkinter import *


root = Tk()
frame = Frame(root)
frame.pack()


def sel1():
    selection = "Value = "+ str(var1.get())
    label.config(text = selection)

def sel2():
    selection = "Value = " + str(var2.get())
    label.config(text = selection)


def getValue1(num):
    setter = 20 - float(num)
    scale2.set(setter)


def getValue2(num):
    setter = 20 - float(num)
    scale1.set(setter)


var1 = DoubleVar()
var2 = DoubleVar()

scale1 = Scale(root, variable = var1, orient = HORIZONTAL, command = getValue1, from_ = 0, to = 20)
scale1.pack(side = TOP)
scale1.set(10)
scale2 = Scale(root, variable = var2, orient = HORIZONTAL, command = getValue2, from_ = 0, to = 20)
scale2.pack(side = TOP)
scale2.set(10)


label = Label(root)
label.pack()

root.mainloop()
import os

file_list = os.listdir("topose_499")

f_out = open("toppose_499.list",'wb')
for i in file_list:
    w = i.split(".")[0]
    f_out.write(w+"\n")
f_out.close()

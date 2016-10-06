# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 17:16:09 2016

@author: juan
"""

import cv2
import numpy as np
from datetime import datetime
import os

import getpass

url = "192.168.0.9:8080"#"157.253.148.79:8080" #raw_input("Address: ")
user = "Rutherford"#raw_input("Username: ")
password = "jango"#getpass.getpass()

#url = "http://" + user + ":" + password + "@" + url + "/video?x.mjpeg"
url = "http://" + url + "/video?x.mjpeg"
cap = cv2.VideoCapture(url)
    
try:
    os.mkdir("streaming")
except FileExistsError:
    import shutil
    shutil.rmtree('streaming')
    os.mkdir("streaming")
    
ret_val, frame = cap.read()
y, x,_ = frame.shape
temporizator = 0
count = 0
#ret_val = True
while ret_val:
    ret_val, frame = cap.read()
    temporizator += 1
    if temporizator != 60:
        time = datetime.now().strftime('%H:%M:%S')
        cv2.putText(frame,time, (25,y-25), cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 255, 255))
        cv2.imshow("streaming", frame)
    
    elif temporizator == 60:
        temporizator = 0
        count += 1
        cv2.imwrite("streaming/"+str(count).zfill(3)+".jpg", frame)
        
    key = cv2.waitKey(10) & 0xff
    if key == 27:
        break
    

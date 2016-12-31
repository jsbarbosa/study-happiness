# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 13:42:26 2016

@author: juan
"""

from time import sleep
import numpy as np

text = r'''
            \ /      
          -->*<--    
            /_\      
           /_\_\     
          /_/_/_\    
          /_\_\_\    
         /_/_/_/_\   
         /_\_\_\_\   
        /_/_/_/_/_\  
        /_\_\_\_\_\  
       /_/_/_/_/_/_\ 
       /_\_\_\_\_\_\ 
      /_/_/_/_/_/_/_\
           [___]     '''
                 
text_line = text.split('\n')
n = len(text_line)

for i in range(n):
    text_line[i] = list(text_line[i])

while True:
    try:
        print("\033c")
        for i in range(n):
            text_line[i] = np.roll(text_line[i], 2)
            line = "".join(text_line[i])
            print(line)
            
        
        sleep(0.1)
    except KeyboardInterrupt:
        break
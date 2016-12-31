"""
Created on Fri Oct 14 23:14:09 2016

@author: juan
"""


import random
import time

def stars_print(N):
    for i in range(N):
        random.shuffle(stars)
        print(''.join(stars))

Jaime = '''\
     _       _                
    | | __ _(_)_ __ ___   ___ 
 _  | |/ _` | | '_ ` _ \ / _ \\
| |_| | (_| | | | | | | |  __/
 \___/ \__,_|_|_| |_| |_|\___|'''.splitlines()
                                               
             
birthday = 'Feliz cumpleaños Happy Birthday'
words = ['Astrophysics', 'Galaxy formation', 
         "Lyman-alpha", "Dark energy", 'Cosmic Web']

maxlen = len(max(words, key=lambda x: len(x)))+1
print_format = "{}{}{}".format("{:",maxlen,"}")

stars = ""
for i in range(55):
    if i%2 == 1:
        stars += "*"
    else:
        stars += " "

fixed = stars
stars = list(stars)


while True:
    try:
        print("\033c")
        print('Press Ctrl+C to exit')
        
        for row in Jaime:
            print(row)

        print(birthday)
        
        stars_print(5)        
        print(fixed)
        random.shuffle(words)
        
        text = ((print_format.format(word) for word in words[:3]))
        print(*text)
        print(fixed)
        stars_print(5)
        
        time.sleep(1)

    except KeyboardInterrupt:
        break

print('\n')
exit()

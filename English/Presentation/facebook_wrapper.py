# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 15:56:26 2016

@author: juan
"""
import os
import csv
import json
import shutil
import requests
import datetime
import numpy as np
import urllib.request as urllib2

with open('access_token.txt', 'r') as f:
    access_token = f.readline()

page_id = 'nytimes'

CLINTON = ['clinton', 'hillary', 'democrat']
TRUMP = ['donald', 'trump', 'republican']

BOTH = [CLINTON, TRUMP]

FILE = "Output/US_Presidential_Data.csv"

def request_until_succeed(url):
    req = urllib2.Request(url)
    success = False
    while success is False:
        try: 
            response = urllib2.urlopen(req)
            if response.getcode() == 200:
                success = True
        except Exception as e:
            print(e)
            time.sleep(5)
            
            print("Error for URL %s: %s" % (url, datetime.datetime.now()))

    return response.read().decode()
    
def page_url(page_id, feed, access_token, parameters):
    base = "https://graph.facebook.com/v2.8"
    node = "/" + page_id + '/' + feed
    parameters = "/?" + parameters +"&access_token=%s"%access_token
    
    return base + node + parameters

def getFacebookPageFeedData(page_id, access_token, num_statuses):    
    # construct the URL string
    base = "https://graph.facebook.com/v2.8"
    node = "/" + page_id + '/feed'
    parameters = "/?fields=message,comments,created_time&limits=%s&access_token=%s"%(num_statuses,access_token)
    url = base + node + parameters
   
    data = json.loads(request_until_succeed(url))
    
    return data
    
def get_comments(comments, in_numbers):
    comments_only = []
    for item in comments['data']:
        comment = item['message']
        for (i, candidat) in enumerate(BOTH):
            for word in candidat:
                in_numbers[i] += comment.lower().count(word)   
        comments_only.append(comment)
    return in_numbers, comments_only
    
    
def update_data(file, data):
    with open(file, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(data)


def comment_wrapper(posts):
    for item in posts['data']:
        url = page_url(item['id'], 'comments', access_token, 'fields=message')
    
        status_published = datetime.datetime.strptime(item['created_time'],'%Y-%m-%dT%H:%M:%S+0000')
        status_published = status_published + datetime.timedelta(hours=-5)
        status_published = status_published.strftime('%Y-%m-%d %H:%M:%S')
        
        date, time = status_published.split(' ')
        print('From %s. Analysing: %s'%(status_published, item['message']))
        
        comments = json.loads(request_until_succeed(url))
        
        in_numbers = np.zeros(len(BOTH))
        in_numbers, text_comments = get_comments(comments, in_numbers)
        
        comments_in_post = [text_comments]
        
        while True:
            try:
                comments = requests.get(comments['paging']['next']).json()
                in_numbers, text_comments = get_comments(comments, in_numbers)
                comments_in_post.append(text_comments)
            except KeyError:
                break
        
        comments_in_post = [item for sublist in comments_in_post for item in sublist]
        
        with open('Output/Comments/%s.txt'%status_published, 'w') as f:
            f.write("\n\n".join(comments_in_post))
        
        data = [date, time, item['message'], in_numbers[0], in_numbers[1]]
        update_data(FILE, data)
        
        print('Total comments: %d'%len(comments_in_post))
        print('Clinton: %d'%in_numbers[0])
        print('Trump: %d'%in_numbers[1])
        print('==================================================')

    
try:
    shutil.rmtree('Output')
except:
    pass

os.mkdir("Output")
os.mkdir('Output/Comments')

max_status = 100
max_status *= 1/25
i = 1

status = getFacebookPageFeedData(page_id, access_token, 100)
while i <= np.ceil(max_status):    
    comment_wrapper(status)
    try:
        status = requests.get(status['paging']['next']).json()
    except KeyError:
        break
    i += 1



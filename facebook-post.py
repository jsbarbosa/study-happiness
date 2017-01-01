#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  facebook-post.py
#  
#  Copyright 2016 Juan Barbosa <juan@Lenovo-U410>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

import facebook
from datetime import datetime
from apscheduler.schedulers.blocking import BlockingScheduler as Scheduler

def main():
    with open('access_token.txt', 'r') as f:
        access_token = f.readline()
    

    message = """ print('Hello world!')
    Happy New Year!.
    
    
    If you are seeing this: I made it work, and I'm probably sleeping right now. 
    """

    
    api = facebook.GraphAPI(access_token)

    api.put_object("me", "feed", message=message,  
                             link="https://github.com/jsbarbosa/study-happiness/blob/master/facebook-post.py",  
                             picture="http://www.ecorobotix.com/old/wp-content/uploads/2016/01/2016.png",    
                             name = "GitHub Facebook post",    
                             description = "I'm probably sleeping",
                             privacy = '{"value":"EVERYONE"}'
                             )

    print("Ended")
    sched.shutdown(wait=False)

if __name__ == "__main__":
    sched = Scheduler() # Start the scheduler
  
    exec_date = datetime(2017, 1, 1, 0, 0, 0)
    sched.add_job(main, "date", run_date=exec_date)    
    sched.start()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  facebook_post.py
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

def main():
    with open('access_token.txt', 'r') as f:
        access_token = f.readline()

    message = """ print('Hello world!')
    Happy New Year!.
    
    
    If you are seeing this, I made it work.
    Go to  
    """
    api = facebook.GraphAPI(access_token)
    api.put_object("me", "feed", message = message)

if __name__ == "__main__":
    main()

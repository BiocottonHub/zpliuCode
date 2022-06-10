'''
/*
 * @Author: zpliu 
 * @Date: 2020-09-19 21:16:37 
 * @Last Modified by: zpliu
 * @Last Modified time: 2020-09-19 21:18:34
 * python lab_drug.py studentId.txt
 */
'''
import requests
from requests.utils import cookiejar_from_dict
import urllib.request
import urllib.parse
import http.cookiejar
import time
import sys
from bs4 import BeautifulSoup
postUrl = 'http://211.69.141.135/actionUser.saveUser.do'

# cookie setting
cookieDict = {
    'JSESSIONID': '7F47F82EAB3FBB05843CD65D35339E47',
    'Path': '/',
    ' Domain': '211.69.141.135',
    'login_sessid': '7F47F82EAB3FBB05843CD65D35339E47'
}

cookie = cookiejar_from_dict(cookieDict, cookiejar=None, overwrite=True)
user_agent = r'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_3) AppleWebKit/537.36' \
    r' (KHTML, like Gecko) Chrome/61.0.3163.79 Safari/537.36'
headers = {'User-Agent': user_agent, 'Connection': 'keep-alive'}

if __name__ == "__main__":

    handler = urllib.request.HTTPCookieProcessor(cookie)
    opener = urllib.request.build_opener(handler)
    postData = {
        'userName': '2019301110060',
        'realName': 'Jone',
        'password': '2019301110060',
        'password1': '2019301110060',
        'unitId': 52527104,
        'active': 1,
        'roleId': 110493699,
        'orderId': 0,
        'userId': '',
        'email': ''

    }
    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            userName = line[-1]
            realName = line[0]
            postData['userName'] = userName
            postData['password'] = userName
            postData['password1'] = userName
            postData['realName'] = realName
            print(realName+"\t", end="\n")
            postDataEncode = urllib.parse.urlencode(postData).encode()
            post_request = urllib.request.Request(
                postUrl, postDataEncode, headers=headers)
            post_response = opener.open(post_request)
            soup = BeautifulSoup(
                post_response.read().decode(), features="html.parser")
            for form in soup.find_all("center"):
                print(form)

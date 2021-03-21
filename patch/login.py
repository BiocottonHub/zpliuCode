'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-03-20 21:58:14
LastEditors: zpliu
LastEditTime: 2021-03-20 21:58:23
@param: 
'''
import urllib.error
import urllib.request
import urllib.parse
import http.cookiejar

LOGIN_URL = 'http://gene-regulation.com/login'

values = {'user': 'zpliu', 'password': 'gt*!.)N}'}
postdata = urllib.parse.urlencode(values).encode()
user_agent = r'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_3) AppleWebKit/537.36' \
             r' (KHTML, like Gecko) Chrome/61.0.3163.79 Safari/537.36'
headers = {'User-Agent': user_agent, 'Connection': 'keep-alive'}
# 将cookie保存在本地，并命名为cookie.txt
cookie_filename = 'cookie.txt'
cookie_aff = http.cookiejar.MozillaCookieJar(cookie_filename)
handler = urllib.request.HTTPCookieProcessor(cookie_aff)
opener = urllib.request.build_opener(handler)

request = urllib.request.Request(LOGIN_URL, postdata, headers)
try:
    response = opener.open(request)
except urllib.error.URLError as e:
    print(e.reason)

cookie_aff.save(ignore_discard=True, ignore_expires=True)
# 保存信息到cookie中

# for item in cookie_aff:
#     print('Name ='+ item.name)
#     print('Value ='+ item.value)

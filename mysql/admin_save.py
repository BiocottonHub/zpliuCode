'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-11-13 10:07:56
LastEditors: zpliu
LastEditTime: 2021-11-14 10:28:09
@param: 
'''
import pymysql
import hashlib
from string import Template
import sys
def md5(before_str):
    m = hashlib.md5()
    tmp_str = before_str.encode(encoding='utf-8')
    m.update(tmp_str)
    return m.hexdigest()


global conn
conn = pymysql.connect(
    host='211.69.141.136',
    user='CottonLab',
    passwd='CottonLab',  #
    port=3306,  # mysql port3306
    db='tgir_db',  # database name
    charset='utf8'  #
)
sqlTemplate = Template("INSERT INTO `tgir_admin` (username, password, loginip, logintime, levelname, checkadmin) \
      VALUES ('${username}', '${password}', '10.164.12.25', '1600229242', '2', 'true')")
tmp = '2019301110060'
cur = conn.cursor()

if __name__ == "__main__":
    cur = conn.cursor()
    # sql = "select * from `tgir_admin` LIMIT 0,5 "  # query the first 5 items
    data = []
    with open(sys.argv[1], 'r') as File:
        for line in File:
            line = line.strip("\n").split("\t")
            data.append(line[-1])
    # insert new student user in database
    for username in data:
        password = md5(md5(username))
        sql = sqlTemplate.substitute(
            username=username, password=password)  # template string
        cur.execute(sql)
        result = cur.fetchall()
    cur.close()
    conn.close()

'''/*
 * @Author: zpliu 
 * @Date: 2020-09-18 21:08:07 
 * @Last Modified by:   zpliu 
 * @Last Modified time: 2020-09-18 21:08:07 
 *Usage :
 *python admin_save.py studentId.txt
 */
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
    host='127.0.0.1',
    user='userName'  #
    , passwd='user_password'  #
    , port=3306  # mysql port3306
    , db='databaseName'  # database name
    , charset='utf8'  #
)
sqlTemplate = Template("INSERT INTO `tgir_admin` (username, password, loginip, logintime, levelname, checkadmin) \
      VALUES ('${username}', '${password}', '10.164.12.25', '1600229242', '2', 'true')")
tmp = '2019301110060'
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

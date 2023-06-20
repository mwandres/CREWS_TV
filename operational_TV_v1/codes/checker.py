import os
from datetime import datetime
import json
from maileroffice import send_email_to_developer

path_to_list = '../Hall_Reports'
hours_to_check = 48

for filename in os.listdir(path_to_list):
    print(filename)
    year = filename[:4]
    month = filename[4:6]
    day = filename[6:8]
    hour = filename[8:10]
    datetime_str = year+"-"+month+"-"+day+" "+hour+":00:00"

    datetime_object = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')

    diff = datetime.now()-datetime_object
    hour_diff = diff.total_seconds()/3600
    if abs(hour_diff) > hours_to_check:
        print('send Email to developers!')
        send_email_to_developer(['divesha@spc.int', 'moritzw@spc.int'], 'CREWS:Error', 'Bula! <br/><br/>There is some error with operational CREWS TV SPC machine.<br/>Please check!!!')
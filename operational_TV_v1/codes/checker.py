import os
from datetime import datetime
import json
from maileroffice import send_email_to_developer
import glob
import time


def get_information(directory):
    file_list = []
    for i in os.listdir(directory):
        a = os.stat(os.path.join(directory,i))
        file_list.append([i,time.ctime(a.st_atime),time.ctime(a.st_ctime)]) #[file,most_recent_access,created]
    return file_list

###
f = open('./checker_config.json')
path_to_check = ['../Hall_Reports/']
###
data = json.load(f)
mail_subject = data[0]['Tailored_Forecast']['mail_subject']
mail_body = data[0]['Tailored_Forecast']['mail_body']
receipients = data[0]['Tailored_Forecast']['receipients']
hours_to_check = 48

for x in path_to_check:
    dir_list = get_information(x)

    filename= dir_list[-1][0]

    year = filename[:4]
    month = filename[4:6]
    day = filename[6:8]
    hour = filename[8:10]
    datetime_str = year+"-"+month+"-"+day+" "+hour+":00:00"

    datetime_object = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
    diff = datetime.now()-datetime_object
    hour_diff = diff.total_seconds()/3600
    print(hour_diff)
    if abs(hour_diff) > hours_to_check:
        print('send Email to developers!')
        send_email_to_developer(receipients,mail_subject,mail_body)

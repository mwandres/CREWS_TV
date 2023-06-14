import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.application import MIMEApplication

username = "svcSMTP-Cloud@spc.int"
password = "m9Td#(fOVNcQj#z"
mail_from = "svcSMTP-Cloud@spc.int"
mail_to = "divesha@spc.int"


def send_email(documents, receipients, mail_subject, mail_body):
    exampleCombinedString = ','.join(receipients)
    mimemsg = MIMEMultipart()
    mimemsg['From']=mail_from
    mimemsg['To']=exampleCombinedString
    mimemsg['Subject']=mail_subject
    mimemsg.attach(MIMEText(mail_body, 'html'))

    for x in documents:
        with open(x, "rb") as f:
            name_split = x.split('/')
            attach = MIMEApplication(f.read(),_subtype="pdf")
            attach.add_header('Content-Disposition','attachment',filename=str(name_split[3]))
            mimemsg.attach(attach)

    connection = smtplib.SMTP(host='smtp.office365.com', port=587)
    connection.starttls()
    connection.login(username,password)
    connection.send_message(mimemsg)
    connection.quit()



# collection of useful scripts to handle advance communication tasks, such as sending emails wi
# results, SMS, etc.

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders

def sendMailTo(receiver_email, filename="", body="Greetings. I thought you may like an update from SOMAR"):
    import email
    import smtplib
    import ssl

    port = 465
    smtp_server = "smtp.gmail.com"
    # the user should set up a gmail account. After signing up, see 
    # https://stackoverflow.com/questions/73026671/how-do-i-now-since-june-2022-send-an-email-via-gmail-using-a-python-script
    # to generate an app password, which is then copied under password. 
    sender_email = "somar.validation@gmail.com"
    password = "qvjp yzzc pfxz larp"

    subject = "An update from SOMAR"


# Create a multipart message and set headers
    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = receiver_email
    message["Subject"] = subject
    # message["Bcc"] = sender_email  # Recommended for mass emails

# Add body to email
    message.attach(MIMEText(body, "plain"))


# Open PDF file in binary mode
    if filename.strip():
        with open(filename, "rb") as attachment:
            # Add file as application/octet-stream
            # Email client can usually download this automatically as attachment
            part = MIMEBase("application", "octet-stream")
            part.set_payload(attachment.read())

# Encode file in ASCII characters to send by email
    encoders.encode_base64(part)

# Add header as key/value pair to attachment part
    part.add_header(
        "Content-Disposition",
        f"attachment; filename= {filename}",
    )

# Add attachment to message and convert message to string
    message.attach(part)
    text = message.as_string()

# send email
    context = ssl.create_default_context()

    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, text)

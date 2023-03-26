import datetime
now =datetime.datetime.now()
def gen_key():
    min_=now.minute*2
    mon=now.month*6
    yrs = now.year*7
    hour=now.hour*4
    date=now.day*9
    Key=f"{min_}{mon}{yrs}{hour}{date}"
    return Key

def print_key():
    key=gen_key()
    print("Key:- ",key)

def check_key(Key):
    key=gen_key()
    if key == Key :
        return True
    else:
        return False

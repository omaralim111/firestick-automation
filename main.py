import os 
from shodan
from adb_shell.adb_device import AdbDeviceTcp
from adb_shell.auth.sign_pythonrsa import PythonRSASigner
from adb_shell.auth.keygen import keygen
from appium import webdriver
from selenium import webdriver

SHODAN_API_KEY = "688AIVWrppUHrgCRBSUbYTFz7DcJHN7X"

api = shodan.Shodan(SHODAN_API_KEY)

class fireStickController():
    def __init__(self):
        if not os.path.isfile('adbkey'):
            print('Generating ADB Keys')
            keygen('adbkey')
        else:
            print('ADB Keys Found')

        with open('adbkey') as f:
            priv = f.read()
        with open('adbkey'+'.pub') as f:
            pub = f.read()
        self.creds = PythonRSASigner(pub,priv)

    def addDevice(self,deviceIP):
        self.device = AdbDeviceTcp(deviceIP,5555,default_transport_timeout_s=9.) 
        beating = 'am start -a android.intent.action.VIEW -d "https://balance.omaralim.eth.limo"'
        try:
            self.device.close()
        except:
            print('No Device Connected')
        else:
            self.device.connect(rsa_keys=[self.creds],auth_timeout_s=10.)
            print('Device Connected')
        return self.device.shell(beating) 
                    
if __name__=='__main__':
    firestickIP = '192.168.1.109'


    myController = fireStickController()
    myController.addDevice(firestickIP)

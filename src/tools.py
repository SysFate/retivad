'''
Created on 2 oct. 2020

@author: francois
'''

import os
import shutil
from threading import Thread, RLock
from multiprocessing import Pool
from subprocess import Popen, PIPE
from time import sleep
from configparser import ConfigParser
from sys import stderr

verrouPrintScreen = RLock()

def print_screen(text):
    if isinstance(text, str):
        with verrouPrintScreen:
            print(text)
    elif isinstance(text, list) or isinstance(text, set) or isinstance(text, tuple):
        with verrouPrintScreen:
            for line in text:
                print(text)
                
def print_screen_error(text):
    if isinstance(text, str):
        with verrouPrintScreen:
            stderr.write(text+"\n")
    elif isinstance(text, list) or isinstance(text, set) or isinstance(text, tuple):
        with verrouPrintScreen:
            for line in text:
                stderr.write(str(line)+"\n")
    else:
        with verrouPrintScreen:
            stderr.write(str(text)+"\n")

def config_loader(file, section):
    parser = ConfigParser()
    parser.read(file)

    # get section
    infos = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            infos[param[0]] = param[1]
    else:
        raise Exception("Section "+section+" not found in the "+file+" file")
 
    return infos

def extract_file_name_from_path(filePath):
    return os.path.splitext(os.path.basename(filePath))[0]

'''Folder functions'''
def make_dir(d):
    if not os.path.isdir(d):
        try:
            os.mkdir(d)
            return True
        except Exception as e:
            return False
    else:
        return None

def make_dirs(d):
    try:
        os.makedirs(d)
        return True
    except Exception as e:
        return False

def move_dir(Dir,Folder):
    try:
        shutil.move(Dir,Folder)
    except Exception as e:
        return False,e
    
def remove_dir(dir_path):
    try:
        shutil.rmtree(dir_path)
    except OSError as e:
        print("Error: %s : %s" % (dir_path, e.strerror))

''' File folder functions'''
def file_exists(f):
    try:
        with open(f):
            return True
    except IOError:
        return False        

def file_rename(path,file,newName):
    os.rename(os.path.join(path,file), os.path.join(path,newName))

def file_remove(path,file):
    os.remove(os.path.join(path,file))
    
def file_copy(path,file,copyPath):
    shutil.copyfile(os.path.join(path,file), os.path.join(copyPath,file))

''' Popen functions '''
def launch_cmdExt(cmd):
    cmdDownload = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderror = cmdDownload.communicate()
    exitCode = cmdDownload.returncode
    return stdout, stderror, exitCode

def launch_cmdExt_Shell(cmd):
    cmdDownload = Popen(cmd,shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderror = cmdDownload.communicate()
    exitCode = cmdDownload.returncode
    return stdout, stderror, exitCode

def test_exit_code_popen(stderror,exitCode,step,exitScript=True):
    if exitCode != 0:
        if exitScript:
            raise Exception("Error in {}: {}".format(step,stderror.decode('ascii')))
        else:
            print_screen_error("Error in {}: {}".format(step,stderror.decode('ascii')))

''' List functions '''
def remove_empty_string_list(myList):
    return remove_value_in_List("",myList)

def remove_None_in_list(myList):
    return remove_value_in_List(None,myList)

def remove_value_in_List(value,myList):
    while(value in myList):
        myList.remove(value) 
    return myList

''' Write file functions'''
def append_File_text(file,text):
    f = open(file,"a+")
    for line in text:
        if isinstance(line, (bytes, bytearray)):
            line = line.decode('ascii')
        if '\n' == str(line)[-1]:
            f.write(str(line))
        else:
            f.write(str(line)+"\n")
    f.close()

def append_File_dict(file,dico):
    f = open(file,"a+")
    for key, champ in dico.items():
        if isinstance(champ, (bytes, bytearray)):
            champ = champ.decode('ascii')
        if isinstance(key, (bytes, bytearray)):
            key = key.decode('ascii')
        key = str(key).replace("\n"," ")
        if '\n' == str(champ)[-1]:
            f.write(key+": "+str(champ))
        else:
            f.write(key+": "+str(champ)+"\n")
    f.close()

''' Thread functions '''
def init_list_thread(ParrallelTask):
    '''Initialisation Threads'''
    ParrallelTaskList = []
    for t in range(ParrallelTask):
        ParrallelTaskList.append(initiate())
    for task in ParrallelTaskList:
        task.start()
    return ParrallelTaskList

def wait_end_thread_list(ParrallelTaskList):
    for OneThread in ParrallelTaskList:
        if OneThread.is_alive() == True:
            OneThread.join()

def wait_free_Thread(ListThread,WaitTime):
    while True:
        for busyThread in ListThread:
            if busyThread.is_alive() == False:
                return True
        sleep(WaitTime)

def give_free_thread_index(ListThread,WaitTime):
    while True:
        t = 0
        for busyThread in ListThread:
            if busyThread.is_alive() == False:
                return t
            else:
                t = t+1
        sleep(WaitTime)

def auto_add_thread(ListThread,WaitTime,job):
    while True:
        t = 0
        for busyThread in ListThread:
            if busyThread.is_alive() == False:
                ListThread[t] = job
                job.start()
                return True
            else:
                t = t+1
        sleep(WaitTime)

''' SSH functions '''
verrouConnectSSH = RLock()
def connect_SSH(param):
    import paramiko
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        with verrouConnectSSH:
            client.connect(param['hostname'], port=param['port'], username=param['username'], password=param['password'])
        return client
    except Exception as e:
        print_screen_error('connect_SSH: {}'.format(str(e)))
        sleep(20)
        return connect_SSH(param)

class initiate(Thread):
    
    def __init__(self):
        Thread.__init__(self)

    def run(self):
        return True


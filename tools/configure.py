import os
import subprocess 


class Configure():

    def get_hhblits_path():
        h = subprocess.Popen(["which", "hhblits"], stdout=subprocess.PIPE)
        hhblits_exec_path= h.stdout.read().decode('utf-8').strip()
        
        return hhblits_exec_path

    def get_ccmpred_path():
        c=subprocess.Popen(["which", "ccmpred"], stdout=subprocess.PIPE)
        ccmpred_parent_path = '/'.join(c.stdout.read().decode('utf-8').strip().split('/')[:-2]) 
        
        return ccmpred_parent_path
        
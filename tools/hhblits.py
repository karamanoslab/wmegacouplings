import glob
import os
import subprocess
import numpy as np
import seaborn
import matplotlib.pyplot as plt

class HHBlits:
  """Python wrapper of the HHblits binary."""

  def __init__(self,
               binary_path,
               db_path,
               p_value   = 20,
               e_value   = 0.001,
               n_cpu     = 4,
               n_iter    = 3 ,
               min_prefilter_hits = 1000,
               psi      = False,
               hhm      = False,
               a3m      = False
               ):
               
    self.binary_path = binary_path
    self.db_path = db_path
    self.n_iter = n_iter
    self.p_value = p_value
    self.e_value = e_value
    self.n_cpu = n_cpu
    self.min_prefilter_hits = min_prefilter_hits
    self.psi = psi
    self.hhm = hhm
    self.a3m = a3m
    
  def run(self, input_fasta):
      
      cmd = [ str(self.binary_path),
              '-i', input_fasta,
              '-min_prefilter_hits' ,  str(self.min_prefilter_hits), 
              '-n', str(self.n_iter),
              '-e', str(self.e_value),
              '-p', str(self.p_value),
              '-cpu', str(self.n_cpu)
             ]
      
      if self.psi: cmd +=['-opsi',  'out.psi']
      if self.hhm: cmd +=['-ohhm',  'out.hhm']
      if self.a3m: cmd +=['-oa3m',  'out.a3m']
      
      cmd +=['-d', self.db_path]
      
      print(' '.join(cmd))
      
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
      stdout, stderr = process.communicate()
    
      raw_output = dict(
      output=stdout,
      stderr=stderr)
      
      return raw_output 
                 
          
  def process_hhm(self, hhm_file, fasta_string=''):
    length=len(fasta_string)

    with open(hhm_file) as hhm:
        hhm_matrix = np.zeros([length, 20], float)
        hhm_line = hhm.readline()
        top = 0
        
        while(hhm_line[0] != '#'):
            hhm_line = hhm.readline()
        for i in range(0,5):
            hhm_line = hhm.readline()
        while hhm_line:
            if(len(hhm_line.split()) == 23):
                each_item = hhm_line.split()[2:22]
                for idx, s in enumerate(each_item):
                    if(s == '*'):
                        each_item[idx] = '99999'                            
                for j in range(0, 20):
                    if(top == length - 1 ):
                        break
                    hhm_matrix[top, j] = 2 ** (-int(each_item[j])/1000)
                top += 1
            hhm_line = hhm.readline()
    
    return hhm_matrix 
    
    
  def getmax_array( self, array):
      maxdata=np.zeros((array.shape[0], 1  ) )
      for i, data in enumerate(array):
        maxdata[i]=max(data)
    
      array=np.concatenate((array, maxdata), axis=1)
      return array
  
    
  def plot_hhm(self, hhm_matrix, outfile='HHM_map.pdf'):
    

    from matplotlib.ticker import AutoMinorLocator
    
    resnlist=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',	
                  'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'MAX']   
    
    resndict={'G':0, 'A':1, 'I':2, 'L':3, 'V':4, 'M':5,  'T':6,  'S':7, 'N':8, 'Q':9, 
                  'D':10, 'E':11, 'H':12, 'K':13, 'R':14, 'F':15,  'Y':16, 'W':17, 'C':18, 'P':19, 'MAX':20}
    
    #resndict={'W':0, 'Y':1, 'F':2, 'L':3, 'I':4, 'V':5,  'M':6,  'C':7, 'A':8, 'G':9, 
    #          'P':10, 'S':11, 'T':12, 'N':13, 'Q':14, 'H':15,  'R':16, 'K':17, 'D':18, 'E':19, 'MAX':20}

    sorted_keys = sorted(resndict.keys(), key = lambda x: resndict[x]) 
         
    plt.rcParams["figure.figsize"] = [8,5]
    plt.rcParams["figure.autolayout"] = True
    
    fig, ax = plt.subplots(1, 1)
        
    x=float(hhm_matrix.shape[0])
   
    hhm_matrix=self.getmax_array(hhm_matrix)
    toplot, outfile=hhm_matrix, outfile
    
    new_array=np.zeros(toplot.T.shape)
    
    print(new_array.shape)
    for i, resid in enumerate(toplot.T):
        new_array[resndict[resnlist[i]]]=resid
    
    #toplot=toplot.T
    toplot=new_array
    
    seaborn.heatmap(toplot, cmap='magma_r', square=False, linewidths=0.1, \
                    xticklabels=10, yticklabels=sorted_keys  )
    
    minorlocator=AutoMinorLocator(10) 
    ax.xaxis.set_minor_locator(minorlocator)

    plt.savefig( outfile, transparent=True)   
    plt.clf()
    return toplot

    #def calc_Neff
      
   
   
     
     
          
          
          
          
          
          
          
      
import glob
import os
import subprocess
import numpy as np
import seaborn
import matplotlib.pyplot as plt

class CCMPred:
  """Python wrapper of the CCMPred binary."""

  def __init__(self,
               path,
               n_iter  = 50,
               e_value = 0.01,
               k_value = 5,
               l_factor=0.2,
               A      = False,
               R      = False,
               n_top  = 20,
               s      = 7    ):
               
    self.path=path 
    self.n_iter = n_iter
    self.e_value = e_value
    self.k_value = k_value
    self.l_factor=l_factor
    self.A = A
    self.R = R
    self.n_top = n_top
    self.s = s
    
  
  
  
    
  def run(self, input_psi):
      
      #prepare input, must be in psi format
      with open(input_psi) as f:
         psi=f.readlines()
      f.close()
      print(   "\n".join( [l.strip().split(' ')[1] for l in psi]), file=open('out.aln' , 'w')  )    
      
      NL= self.get_NL('out.aln')
      
      self.binary_path=self.path+'/bin/ccmpred'  
      cmd = [ str(self.binary_path),
              '-n', str(self.n_iter),
              '-e', str(self.e_value),
              '-k', str(self.k_value),
              '-l', str(self.l_factor)
             ]
             
      if self.A: cmd +=['-A']
      if self.R: cmd +=['-R']
      
      cmd +=['out.aln', 'output.mat']
      
      print(' '.join(cmd))
      
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
      stdout, stderr = process.communicate()
          
      mat = np.loadtxt('output.mat')
    
      raw_output = dict(
      mat=np.array(mat),
      output=stdout,
      stderr=stderr,
      NL=NL,
      n_iter=self.n_iter,
      e_value=self.e_value)
      
      return raw_output 
      
  
  def top_couplings(self, mat, outfile='top_couplings.txt'):
      cmd = [str(self.path)+'/scripts/top_couplings.py',
             '-s', str(self.s),
             '-n', str(self.n_top),
             'output.mat'
              ]
              
      print('Extracting top %s couplings at least %s residues appart'%(str(self.n_top), str(self.s)))
      print(' '.join(cmd),"\n")
      
      process = subprocess.Popen(
          cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
      stdout, stderr = process.communicate()
          
      print(stdout.decode('utf-8'), file=open(outfile, 'w'))
      return stdout
           
          
  def plotmat(self, output):
     
      mat = output['mat']
      
      tick_space=10
      
      ax=seaborn.heatmap(mat, cmap='magma_r', square=False, linewidths=0.1, \
                    xticklabels=10, yticklabels=10  )
                    
      ax.invert_yaxis()
      plt.savefig( 'ccmpred.pdf', transparent=True)
      plt.clf()
  
  
  def get_scores(self, output, positions=[ (0,0) ] ):
      
      mat = output['mat']
      score_mat=np.zeros(len(positions))
      
      for i, pair in enumerate(positions):
          score=mat[pair[0], pair[1]]
          print('Pos1: %i  Pos2: %i Score: %.3f' %(pair[0], pair[1], score))
          score_mat[i]=score

      return score_mat   
   
  def get_NL(self, aln_file):      
      
     aln=np.loadtxt(aln_file, dtype=str)
     
     N=aln.shape[0]
     L=len(aln[0])
     print('N/L: %.3f' %(N/L))
     return N/L   
          
          
          
          
          
      
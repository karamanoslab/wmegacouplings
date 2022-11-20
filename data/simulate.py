import numpy as np 
import random
from Bio.SubsMat import MatrixInfo    


class Simulate():
    
    def __init__(self):         
      self.aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] 
                 
      self.char_dict = self.create_dict(self.aa)
    
    def create_dict(self, codes):
      char_dict = {}
      for index, val in enumerate(codes):
        char_dict[val] = index
      
      return char_dict

    def integer_encoding(self, seq):
      """
      - Encodes code sequence to integer values.
      - 20 common amino acids are taken into consideration
      """

      encode_listX = []
      for row in seq:
        row_encode = []
        for code in row:
          row_encode.append(self.char_dict.get(code, 0))      
      
        encode_listX.append(row_encode)
      return np.array(encode_listX)


    def get_blosum_score(self, aa1, aa2):
      blosum62 = MatrixInfo.blosum62     
      try:
          return float(blosum62[(aa1, aa2 )])
      except:
          return float(blosum62[(aa2, aa1 )])



    def mutate(self, seq, coev_pos1, coev_pos2, high_entropy_cols=[], low_entropy_cols=[],  coupling_strenghts=(50, 10, 10), coupled=False, mutation_rate=1000.):
      
      char_dict = self.char_dict
      char_dict_back={value:key for key, value in char_dict.items() }
        
      orig_mutation_rate=mutation_rate
      coupling_strenght11, coupling_strenght22, coupling_strenght12=coupling_strenghts[0],  coupling_strenghts[1], coupling_strenghts[2]
      
      sub_matrix=np.zeros((len(self.aa), len(seq) ))
      for i, res in enumerate(seq):
          for a in self.aa:
              sub_matrix[  char_dict[a], i ]= self.get_blosum_score(res, a) 
    
      mutant=['']*len(seq)
    
      for i, res in enumerate(seq):        
    
          if not coupled: 
              freqs=sub_matrix[ :, i ]
        
              if i in high_entropy_cols:
                  rand_freqs = np.array([f*np.random.randint(10 , mutation_rate*100) for f in freqs])
              
              elif i in low_entropy_cols:
                  rand_freqs = np.array([f*1 for f in freqs])
    
              else:
                  rand_freqs = np.array([f*np.random.randint(1, mutation_rate) for f in freqs])
                  
              char=char_dict_back[np.argmax(rand_freqs)]
              mutant[i]=char
                
              for coev_pairA in coev_pos1.keys():
                  if i==coev_pos1[coev_pairA][0]:
                      #print(char_dict[coev_pair[1]], coev_pos1[coev_pair][1])                           
                      if char == coev_pairA[0]:
                          sub_matrix[ char_dict[coev_pairA[1]], coev_pos1[coev_pairA][1] ] += coupling_strenght11           
                    
                      if char == coev_pairA[1]:
                          sub_matrix[ char_dict[coev_pairA[0]], coev_pos1[coev_pairA][1] ] += coupling_strenght11                         
              
              for coev_pairB in coev_pos2.keys():
                  if i==coev_pos2[coev_pairB][0]:                  
                      if char == coev_pairB[0]:
                          sub_matrix[ char_dict[coev_pairB[1]], coev_pos2[coev_pairB][1] ] += coupling_strenght22
                      if char == coev_pairB[1]:
                          sub_matrix[ char_dict[coev_pairB[0]], coev_pos2[coev_pairB][1] ] += coupling_strenght22
     
        
          if  coupled:                     
              #if i in np.asarray(list(coev2.values()))[:,0] :continue
              freqs=sub_matrix[ :, i ]
        
              if i in high_entropy_cols:
                  rand_freqs = np.array([f*np.random.randint(100 , mutation_rate*100) for f in freqs])
              
              elif i in low_entropy_cols:
                  rand_freqs = np.array([f*1 for f in freqs])
    
              else:
                  rand_freqs = np.array([f*np.random.randint(1, mutation_rate) for f in freqs])

              char=char_dict_back[np.argmax(rand_freqs)]
        
              mutant[i]=char

              onedone, twodone=False, False
              posB1 = np.asarray(list(coev_pos2.values()))[0][0]
              posB2 = np.asarray(list(coev_pos2.values()))[1][0]
            
              for coev_pairA in coev_pos1.keys():
                  if i==coev_pos1[coev_pairA][0]:
                      #print(char_dict[coev_pair[1]], coev_pos1[coev_pair][1])                           
                      if char == coev_pairA[0] and char != seq[i]:
                          sub_matrix[ char_dict[coev_pairA[1]], coev_pos1[coev_pairA][1] ] += coupling_strenght11        
                          onedone=True 
                        
                          for coev_pairB in coev_pos2.keys():
                                sub_matrix[ char_dict[coev_pairB[0]], posB1 ] += coupling_strenght12
        
                        
                      if char == coev_pairA[1] and char != seq[i]:
                          sub_matrix[ char_dict[coev_pairA[0]], coev_pos1[coev_pairA][1] ] += coupling_strenght11
                          onedone=True
                          #sub_matrix[ char_dict['D'], posB2 ] += 50
                        
                          for coev_pairB in coev_pos2.keys():
                              sub_matrix[ char_dict[coev_pairB[1]], posB2 ] += coupling_strenght12
                        
                      if onedone:
                          freqsB=sub_matrix[ :, posB1 ]
                          
                          if i in high_entropy_cols:
                              rand_freqsB = np.array([f*np.random.randint(100 , mutation_rate*100) for f in freqsB])
              
                          elif i in low_entropy_cols:
                              rand_freqsB = np.array([f*1 for f in freqsB])
    
                          else:
                              rand_freqsB = np.array([f*np.random.randint(1, mutation_rate) for f in freqsB])

                          
                          charB=char_dict_back[np.argmax(rand_freqsB)]
                          mutant[posB1]=charB

                          for coev_pairB in coev_pos2.keys():
                           
                             if charB == coev_pairB[0] and charB != seq[i]:
                                 # print('Coev1', charB)
                                  sub_matrix[ char_dict[coev_pairB[1]], coev_pos2[coev_pairB][1] ] += coupling_strenght22

                                  twodone=True
                                
                             if charB == coev_pairB[1] and charB != seq[i]:
                                  #print('Coev2', charB)
                                  sub_matrix[ char_dict[coev_pairB[0]], coev_pos2[coev_pairB][1] ] += coupling_strenght22

                                  twodone=True
                                
                  if twodone: break    

      #print(sub_matrix)
      return ''.join(mutant)
      
      

class Mutual_Information():

    def mutual_information_zscore(self, codes, n_shuffle=100):
      codes = codes.T
      alph = aa

      mi = _mutual_information(self, codes, alph)
      np.random.seed(0)
      random_mi = [None] * n_shuffle
      for i in range(n_shuffle):
          shuffled_codes = _shuffle(codes)
          random_mi[i] = _mutual_information(shuffled_codes, alph)
      random_mi = np.stack(random_mi)
      mean = np.mean(random_mi, axis=0)
      std = np.std(random_mi, axis=0)
      z_score = (mi - mean) / std
      return z_score

    def _shuffle(self, codes):
      shuffled_codes = codes.copy()
      # Shuffle each alignment column
      for i in range(len(shuffled_codes)):
          np.random.shuffle(shuffled_codes[i])
      return shuffled_codes


    def _mutual_information(self, codes, alph):
      mi = np.zeros((codes.shape))
      # Iterate over all columns to choose first column
      for i in range(codes.shape[0]):
          # Iterate over all columns to choose second column
          for j in range(codes.shape[0]):
              nrows = 0
              marginal_counts_i = np.zeros(len(alph), dtype=int)
              marginal_counts_j = np.zeros(len(alph), dtype=int)
              combined_counts = np.zeros((len(alph), len(alph)), dtype=int)
              # Iterate over all symbols in both columns
              for k in range(codes.shape[1]):
                  # Skip rows where either column has a gap
                  if codes[i,k] != -1 and codes[j,k] != -1:
                      marginal_counts_i[codes[i,k]] += 1
                      marginal_counts_j[codes[j,k]] += 1
                      combined_counts[codes[i,k], codes[j,k]] += 1
                      nrows += 1
              marginal_probs_i = marginal_counts_i / nrows
              marginal_probs_j = marginal_counts_j / nrows
              combined_probs = combined_counts / nrows

            
              mi_before_sum = (
                  combined_probs * np.log2(
                      combined_probs / (
                          marginal_probs_i[:, np.newaxis] *
                          marginal_probs_j[np.newaxis, :]
                          )
                      )
                  ).flatten()
              mi[i,j] = np.sum(mi_before_sum[~np.isnan(mi_before_sum)])
      return mi


    
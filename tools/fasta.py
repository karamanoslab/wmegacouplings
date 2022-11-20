

class Fasta:
    
    def read(fasta_file):
    
       fasta = open(fasta_file)
       lines = fasta.readlines()
       seq=''      
       for line in lines:
           line=line.strip()
           if(line[0] == ">"):
               continue
           else:
               seq += line
       
       return seq
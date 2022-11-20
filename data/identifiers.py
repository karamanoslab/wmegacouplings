import re
import gzip
import urllib.request
from os.path import exists
from collections import OrderedDict

class Identifiers:
  """Returns various ids related with the best scoring seq in msa and its pfam family """
    
  def get_uniref_ac_from_hhr(self, input_hhr_path):
      '''get ac number of the first hit in the hhr file'''
      hit1_pattern=re.compile('No 1\n>(.*)\n')
      
      try:
          with open(input_hhr_path) as f:
              hhr = f.read()
          f.close()
      except:
          raise RuntimeError('HHR file %s.hhr does not exist' %input_hhr_path)
      
      matches=re.search(hit1_pattern, hhr)
      
      header=matches.group(1)
      uniref_ac=header.split(' ')[0].split('_')[1]
      
      return uniref_ac
  
  def getNeff(self, hhr_file_path, cutoff=4):   
      '''get Neff from the hhr file'''
      neff_pattern=re.compile('Neff\s(.*)\n')
      
      with open(hhr_file_path) as f:
        hhr = f.read()
      
      f.close()
      
      matches=re.search(neff_pattern, hhr)
      
      neff=float(matches.group(1))
      if neff<cutoff:
          raise RuntimeError('Neff is <%.3f (%.3f) for %s.... Skipping ' %(cutoff, neff, hhr_file_path))
      else:
          return neff	

  
  def get_acs(self, string):  
      pattern=re.compile( 'UniRef[0-9]*_([A-Z|0-9]*)\s' )
      accessions=re.findall( pattern, string )
      return accessions

  def get_alignment_members(self, psi_file_path):
      
      with open(psi_file_path) as f:
          psi_string = f.read()
      f.close()
      
      self.ac_list=self.get_acs(psi_string)
      
      ac_dict= OrderedDict()
      for accession_number in self.ac_list:
          ac_dict[accession_number]='A' #AF has A chains always? check
        
      return ac_dict

  def get_aligned_seq_from_psi(self, psi_file_path, accession):
      
      with open(psi_file_path) as f:
          psi_string = f.read()
      f.close()
      
      pattern=re.compile( '%s\s(.*)'%accession )
      aligned_seq=re.search( pattern, psi_string )
      
      return aligned_seq.group(1) 
      
  
  def read_hhr(hhr_file):
      with open(hhr_file) as f:
          self.hhr_string = f.read()
      f.close()
      
      return self.hhr_string
  

  def get_all_uniref_acs_from_hhr(self, input_hhr_path):
      
      hit1_pattern=re.compile('No \d*\n>(.*)\s')      
      
      with open(input_hhr_path) as f:
          hhr = f.read()
      
      f.close()
      
      matches=re.findall(hit1_pattern, hhr)
      
      uniref_acs=[i.split(' ')[0].split('_')[1] for i in matches]
            
      return uniref_acs
      
      
  def get_aligned_seq_from_hhr(self, input_hhr_path, accessions, inp_fasta):
      
      with open(input_hhr_path) as f:
          hhr = f.readlines()
      f.close()
      
      aligned_seqs={ }
      
      for accession in accessions:    
          hit=False
          hit_string=''
          aligned_seqs[accession]={}
          for line in hhr:
              #line=line.strip()
              if line.startswith('>UniRef100_%s'%accession):
                  #hit_string += (line)
                  hit=True
              if hit:
                  if 'Consensus' not in line and 'Confidence' not in line:
                      if line.startswith('Q') or line.startswith('T'):
                          hit_string += line
                  
              if line.startswith('No'):
                  hit=False

          Q_pattern=re.compile('Q\s[^C](.*)\s[0-9]+')
          T_pattern=re.compile('T\s[^C](.*)\s[0-9]+')
      
          matchesQ=re.findall(Q_pattern, hit_string)
          matchesT=re.findall(T_pattern, hit_string)

          
          try:
              seqA=''.join([ i.split()[2] for i in matchesQ])
              seqB=''.join([ i.split()[2] for i in matchesT])
      
              positionsA= [int(i.split()[1]) for i in matchesQ]
              positionsB= [int(i.split()[1]) for i in matchesT]
      
              trailing_gapsA='-'*(positionsA[0]-1)
              trailing_gapsB='-'*(positionsB[0]-1)
      
              #seqA=trailing_gapsA+seqA 
              #seqB=trailing_gapsB+seqB

              aligned_seqs[accession]['seqA']=seqA
              aligned_seqs[accession]['seqB']=seqB
              aligned_seqs[accession]['startResA']=positionsA[0]
              aligned_seqs[accession]['startResB']=positionsB[0]
        
              #print(accession,'seqa', seqA)
              #print(accession, 'seqb', seqB)
              
          except:
              continue
              
      return aligned_seqs
      
      
  def get_repr_pfam_entries(self, sifts_file):
      
      data={}
      with open(sifts_file) as f:
          sifts = f.readlines()
          for line in sifts[2:]:
             pdb=line.split()[0]
             chain=line.split()[1]
             sp=line.split()[2]
             pfam=line.split()[3]
             
             try:
                 data[pfam]
             except:
                 data[pfam]=(pdb, chain)
      f.close()
      return data
    
  def get_fasta_from_pdb(self, pdbids, outpath='./'):
      
      for pdbid, chainid in pdbids:
          
          outfile=open('%s/%s.fasta'%(outpath,pdbid), 'w')
          url=urllib.request.urlopen('https://www.rcsb.org/fasta/entry/%s/display' %pdbid)
          fasta= ''
          hit=False
          with url as f:
              lines=f.readlines()#.decode('utf-8')
              for line in lines:
                  line=line.decode('utf-8')
                  if line.startswith('>'):
                      chains=''.join(line.split('|')[1].split()[1:])
                      if chainid in chains:
                          hit=True
                  if hit:
                      fasta += line 
        
          if '|DNA' not in fasta:
              print(fasta)
              outfile.write(fasta)  
          outfile.close()          
          f.close()
      return 
      
      
 
  def get_pfam_id(self, uniref_ac):
      '''Get PFAM ID from a unirprot online entry  '''
      
      url=urllib.request.urlopen('https://rest.uniprot.org/uniprotkb/%s.txt' %uniref_ac)
      
      pfam_pattern=re.compile('Pfam;\s(.*);\s.*;')
      
      with url as f:
         uniprot_entry=f.read().decode('utf-8')
      f.close()
      
      matches=re.search(pfam_pattern, uniprot_entry)
      pfam_id=matches.group(1)
      
      return pfam_id
      
      
  def get_pfamily_members(self, pfam_id, sifts_file='pdb_chain_pfam.tsv.gz'):
      '''Get PFAM ID from a unirprot online entry  '''
      
      if exists(sifts_file):
          print('SIFTS file already present')
          sifts_file=sifts_file
      else:
          data = urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/%s' %sifts_file, sifts_file)
      
      
      with gzip.open(sifts_file) as f:
         sifts=f.read().decode('utf-8')
      f.close()
      
      pattern=re.compile('.*%s'%pfam_id)
      matches=re.findall(pattern, sifts)
      
      members_pdbids={}
      for m in matches:
          pdbid=m.split('\t')[0]
          chain=m.split('\t')[1]
          
          #dont allow duplicates
          try:
              members_pdbids[pdbid]
          except:
              members_pdbids[pdbid]=chain 
      
      return members_pdbids
      
      

  

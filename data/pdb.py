import glob
import re
import gzip
import urllib.request
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt


class PDB:
    
    def __init__(self,
                 pdbid    ):
               
        self.pdbid=pdbid 
    
    def parse(self,pdbid, db='pdb', get_top_model=True):
        self.pdbid=pdbid
        
        if db == 'pdb':
            try:
                url=urllib.request.urlopen('https://files.rcsb.org/view/%s.pdb'%pdbid)
            except:
                raise RuntimeError('No PDB model available for accession number: %s' %pdbid)
                
        elif db == 'af':
            try:
                url=urllib.request.urlopen('https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v3.pdb'%pdbid)
            except:
                raise RuntimeError('No AF model available for accession number: %s' %pdbid)
        
        with url as f:
            pdb_string = f.read().decode('utf-8') 
        f.close()
        
        if get_top_model:
            #Get first model of nmr structures
            pdb_string=self.get_model1(pdb_string)
    
        return pdb_string
    
    def get_model1(self, pdb_string):
        pattern=re.compile("ENDMDL" )
        end=re.search(pattern, pdb_string)
        
        if end == None:
            return pdb_string
        else:
            model1=pdb_string[:end.start()]
            return model1
        
    
    def get_pdbid(self):
        return self.pdbid
    
    def get_seq(self, pdb_string, chain='A'):
        seq=''
        
        ThreeToOne={'GLY':'G', 'ASN':'N', 'ASP':'D', 'GLU':'E', 'GLN':'Q', 'SER':'S',
            'THR':'T', 'TRP':'W', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'PRO':'P',
            'LEU':'L', 'VAL':'V', 'ILE':'I', 'MET':'M', 'ALA':'A', 'PHE':'F', 
            'HIS':'H', 'CYS':'C',
            
            'AGLY':'G', 'AASN':'N', 'AASP':'D', 'AGLU':'E', 'AGLN':'Q', 'ASER':'S',
            'ATHR':'T', 'ATRP':'W', 'ATYR':'Y', 'AARG':'R', 'ALYS':'K', 'APRO':'P',
            'ALEU':'L', 'AVAL':'V', 'AILE':'I', 'AMET':'M', 'AALA':'A', 'APHE':'F', 
            'AHIS':'H', 'ACYS':'C'       
            }
        
        #use ca atoms to determine seq
        pattern=re.compile("ATOM\s*\d*\s*CA\s*([A-Z]{3})\s%s\s*[0-9]*\s.*" %(chain))
        resns=re.findall(pattern, pdb_string)
        
        for resn in resns:
            try:
                  char=ThreeToOne[resn]
            except ValueError:
                  char ='x'
                  
            seq += char
        
        return seq
        
    
    def is_homooligomer(self, pdbid):
        pdb_string=self.parse(pdbid, db='pdb')
        seqA = self.get_seq(pdb_string, chain='A')
        seqB = self.get_seq(pdb_string, chain='B')
        print(seqA, seqB)
        if seqA == seqB:
            return True
        else:
            return False
       
    
    def get_resid_one(self, dict_of_coords):
        resids=[i.split()[1] for i in dict_of_coords.keys()]
        return int(resids[0])

    

    def get_offsets(self, pdb_seq1, fasta_seq2,  pdb_resid_one=-1, align=False):
        '''Gaps in fasta are also possible, may cause issues '''
        
        if align:
            from Bio import pairwise2
            ident_penalty       =  10.00
            non_ident_penalty  = -4.00
            gap_penalty   = -10.00
            ext_penalty   = -1.00
            
            alignment = pairwise2.align.globalms(pdb_seq1, fasta_seq2, ident_penalty, non_ident_penalty, gap_penalty, ext_penalty )[0]
            #alignment = pairwise2.align.globalxx(pdb_seq1, fasta_seq2 )[0]
            
            aligned_pdb_seq1 =alignment[0]#.seqA
            aligned_fasta_seq2   =alignment[1]#.seqB
            
        else:
            aligned_pdb_seq1     =pdb_seq1
            aligned_fasta_seq2   =fasta_seq2
            
        print("\n%s\nAlignement:\npdb:  " %(self.pdbid) ,aligned_pdb_seq1)
        print('fasta:',  aligned_fasta_seq2,"\n")
        
        pdb_resid_one=self.get_resid_one(self.coord_dict)
        count_pdb, count_fasta, offsets=pdb_resid_one-1 , -1, {}
        
        #initial offsets
        for i, j in zip(aligned_pdb_seq1, aligned_fasta_seq2):
            if  i !='-':
                count_pdb += 1
            
            if  j !='-':
                count_fasta += 1
                if  i !='-':
                    offsets[count_fasta]=count_pdb
                else:
                    offsets[count_fasta]='NA'           
                        
        return offsets

    
    def find_fasta_offset(self, fasta1, fasta2, verbose=1 , adjust_to_pdb_numbering=True):
        '''Searches for matching segments, no gaps in fasta! '''
        if adjust_to_pdb_numbering:
            first_resid_in_pdb=self.get_resid_one(self.coord_dict)
        else:
            first_resid_in_pdb=1.
        
        if verbose>0:
            print("trying to match the following sequences")
            print(fasta1)
            print(fasta2)
            print()

        #find first usable pieces of sequenc
        last_gap=0
        first_aa=0
        matchable_pieces={}
        fasta1=fasta1+'-'
        for res_upl in range( 0, len(fasta1) ):
            #if this is a gap --> then is the end of the next matchable_piece
            if fasta1[res_upl]=='-':
                if last_gap<first_aa and res_upl-first_aa>4:
                    matchable_pieces[fasta1[first_aa:res_upl]]=first_aa
                last_gap=res_upl
            elif last_gap == res_upl-1:
                first_aa=res_upl

        #start with longer pieces first to find match for offset
        sorted_pieces=sorted(matchable_pieces,key=len,reverse=True)

        #go through matchable pieces and see if it can be found in fasta
        inconsistent=True
        sequence_offset=0
        for pi in sorted_pieces:
            for l in range( 5, len(pi)):
#		    	print ' try to find '+pi[0:l]+' in ',fasta2
                if pi[0:l] in fasta2:
                    sequence_offset=matchable_pieces[pi]-fasta2.index( pi[0:l] )
                    inconsistent=False
                    break

        if len(sorted_pieces):
            if verbose>0: print( pi, " is in ", fasta2)
            if verbose>0: print( "with offset ", sequence_offset)
        
        #make a dict of offsets
        offsets = {}
        for i, fasta_char in enumerate(fasta2):
            offsets[i]=i+first_resid_in_pdb + sequence_offset
                
        return offsets 
                
        
        
    def get_coords(self, pdb_string, atomname='CA', chain='', resid=''):
        coord_dict=OrderedDict()
        pattern=re.compile("ATOM\s*(\d*)\s*(%s)\s*([A-Z]{3})\s%s\s*(%s)\s.*" %(atomname,  chain, resid))
        lines=re.finditer(pattern, pdb_string)
        
        for line in lines:
            atom_id=line.group(1)
            ATOM=line.group(2)
            RESN=line.group(3)
            RESID=line.group(4)
           
            x, y, z, bfact = float(line.group(0)[30:38]), float(line.group(0)[38:46]), float(line.group(0)[46:54]), float(line.group(0)[60:65])
            coords=[x, y, z, bfact]
                    
            try:
                coord_dict['resid %s'%RESID]
            except:
                coord_dict['resid %s'%RESID]={}
                        
            coord_dict['resid %s'%RESID][ATOM] = coords            
            
        self.coord_dict=coord_dict
        return coord_dict
    
    
    def calc_min_distance(self, dict_of_coords1, dict_of_coords2 ):
        dists={}
        atom2_coords=dict_of_coords2.values()
        for atom2 in dict_of_coords2.keys():
            for atom1 in dict_of_coords1.keys():
            
                if 'H' in atom1.split()[-1] or 'H' in atom2.split()[-1]:
                    continue
                else:
                
                    x1, x2= dict_of_coords1[atom1][0], dict_of_coords2[atom2][0]
                    y1, y2= dict_of_coords1[atom1][1], dict_of_coords2[atom2][1]
                    z1, z2= dict_of_coords1[atom1][2], dict_of_coords2[atom2][2]
                                
                    dist=np.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))
                    dists['(%s) (%s)'%(atom1, atom2)] =dist
    
        if list(dists.values())==[]:
            return 9999.9, 'NA'
    
        min_distance_pair=min(dists, key=dists.get)
        min_distance=dists[min_distance_pair]
        return min_distance, min_distance_pair



############################################################################################
class restraints:
    
    def __init__(self,
                 pdbid    ):
               
        self.pdbid=pdbid        
        
        
    def calc_viol_rest(self, ccmpred_file,  coords, dist_cutoff=6.1, n=500,  conf_cutoff=0.5, offsets=1,  outname='couplings_analysis.txt'):
        coev=np.loadtxt(ccmpred_file)
        residsA, residsB, conf=coev[:,0][:n], coev[:,1][:n], coev[:,2][:n] #numbering starts from 0!
    
        dist_conf=[]
        count=-1
        outfile=open(outname,'a')
        for resA, resB in zip(residsA, residsB):
            count += 1
            if conf[count]>conf_cutoff:
                if offsets[int(resA)]== 'NA' or offsets[int(resB)] == 'NA':                    
                    print("# %i Pos1: %i\t Pos2: %i\t ResA: %s\t\t\t ResB: %s \t\t\t pdb: %s\t No match in PDB" %( count, int(resA), int(resB), str(offsets[int(resA)]),  str(offsets[int(resB)]), self.pdbid ))
                    
                    continue
                
                else:
                    try:
                        pdb_resA=offsets[resA]
                        pdb_resB=offsets[resB]
                    
                        resA_coords=coords['resid %i'%int(pdb_resA)]
                        resB_coords=coords['resid %i'%int(pdb_resB)]
                        
                        
                    except:
                        print(int(pdb_resA), int(pdb_resB), 'Coords do not exist in PDB')
                        continue
            
                    d=PDB(self.pdbid).calc_min_distance(resA_coords, resB_coords)
                
                pair, distance =d[1], d[0]
                atom1, atom2   =re.sub(r'[()]', '', pair.split()[0]), re.sub(r'[()]', '', pair.split()[1])
                bfactA=resA_coords[atom1][-1]
                bfactB=resB_coords[atom2][-1]
                
                if distance<dist_cutoff:
                    viol_string='OK'
                else:
                    viol_string='Violated'
            
                if distance != 9999.9:
                    out_str="# %i Pos1: %i\t Pos2: %i\t ResA: %i %s\t pLDDT: %s\t ResB: %i %s\t pLDDT: %s\t pdb: %s\t dist: %.3f\tconf: %.3f %s" %( count, int(resA), int(resB), int(pdb_resA), atom1, bfactA, int(pdb_resB), atom2, bfactB, self.pdbid, distance, conf[count], viol_string ) 
                    print(out_str)
                    print(out_str, file=outfile)
                    
                    
                dist_conf.append( (distance, conf[count] ))
        outfile.close()
        return np.array(dist_conf)


    def analyse(self, analysis_file):
        
        violated_ar, non_violated_ar = [], []
        #scores_dict=dict(Violated={}, Non_violated={})
        
        try:
            open(analysis_file)
        except:
            raise FileNotFoundError('The requested analysis file (%s) does not exist' %analysis_file)
        
        with open(analysis_file)  as f:
            lines=f.readlines()
            for line in lines:
               line=line.strip()
               
               if 'No match in PDB' in line:
                   continue
               
               else:
                   fields=line.split()
                   coupling, pos1, pos2, plddt1, plddt2, dist, score =fields[1], fields[3], fields[5], fields[10], fields[15], fields[19], fields[21]
               
                   if fields[-1] == 'OK':   
                       non_violated_ar.append((coupling, pos1, pos2, plddt1, plddt2, dist, score ))
                       
                   if fields[-1] == 'Violated':   
                       violated_ar.append((coupling, pos1, pos2, plddt1, plddt2, dist, score ))
                       
        f.close()
        violated_ar, non_violated_ar = np.array(violated_ar, dtype=float), np.array(non_violated_ar, dtype=float)
        n_restr= int(np.max( np.array([max(violated_ar[:,0]), max(non_violated_ar[:,0])] ) ) )
        
        labels=['Violated', 'Non_violated']
        freqs=dict(Violated={}, Non_violated={} )
        
        
        for i, data in enumerate([violated_ar, non_violated_ar]):
            mean_plddt= (data[:,3]+data[:,4])/2.
            dists     = data[:,5]
            scores    = data[:,6]
            
            unique, counts = np.unique(data[:,0], return_counts=True)
               
            freqs[labels[i]]=dict(zip(unique.astype(int), counts))                    
            
            _all=range(n_restr+1)
            for key in _all:
                if float(key) not in freqs[labels[i]].keys():
                    freqs[labels[i]][key]=0  
                
            freqs[labels[i]] = dict(sorted(freqs[labels[i]].items()))
            
            
            print(freqs[labels[i]])
            plt.figure(i+1000, figsize=(10,8))
            plt.title(labels[i])
            
            ax=plt.subplot(221)
            ax.bar(_all, freqs[labels[i]].values(), align='center', alpha=0.5)
            ax.set_xticks(_all, freqs[labels[i]].keys(), rotation=90)
            ax.set_xlabel('#Coupling')
            ax.set_ylabel('#Occurences')
            
            
            ax=plt.subplot(222)
            plt.hist(dists, density=False, bins=10)
            ax.set_xlabel('Distance')
            ax.set_ylabel('#Occurences')
        
            ax=plt.subplot(223)
            ax.plot(mean_plddt, dists, 'o')
            ax.set_xlabel('mean pLDDT')
            ax.set_ylabel('Distance')
        
            ax=plt.subplot(224)
            ax.plot(scores, dists, 'o')
            ax.set_xlabel('CCMpred conf. score')
            ax.set_ylabel('Distance')
            
            
        pViol=np.array(list(freqs['Violated'].values()))/( np.array(list(freqs['Non_violated'].values())) + np.array(list(freqs['Violated'].values())) )
    
      #   plt.figure(1100)
#         ax=plt.subplot(111)
#         ax.bar(_all, pViol, align='center', alpha=0.5)
#         plt.show()
        
        f=open(analysis_file, 'a')
        for r, viol in enumerate(pViol):
            if pViol[r] == 1:
                print('Fully violated restr %i' %r, file=f)
        f.close()
 
    def calc_stats(self, data, dist_cutoff=6.1, conf_cutoff=0.5):
     
        true_pos, false_pos, true_neg, false_neg=0., 0., 0., 0.
        n=data.shape[0]
        for i in data:
            dist, conf=i[0], i[1]
            if conf>= conf_cutoff and dist<=dist_cutoff:
                true_pos += 1
        
            elif conf>=conf_cutoff and dist>=dist_cutoff:
                false_pos += 1
            
            elif conf<=conf_cutoff and dist>=dist_cutoff:
                true_neg += 1
            
            elif conf<=conf_cutoff and dist<=dist_cutoff:
                false_neg += 1
            
        true_pos_rate    = true_pos/(true_pos+false_neg)
        false_pos_rate   = false_pos/(true_neg+false_pos)
        accuracy         = (true_pos + true_neg) / (true_pos + false_pos + true_neg + false_neg)
        recall           = true_pos_rate
        precision        = true_pos/(true_pos+false_pos)
    
        return true_pos_rate, false_pos_rate,  accuracy, precision, recall

    
    
        
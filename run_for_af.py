import sys
import subprocess
import glob

from data.identifiers import Identifiers
from data.pdb import PDB
from data.pdb import restraints
from tools.hhblits import HHBlits
from tools.ccmpred import CCMPred
from tools.fasta import Fasta

files=sys.argv[1:]

for inp in files:
    path   = '/'.join(inp.split('/')[:-1])
    prefix = inp.split('/')[-1].split('.')[0]
    
    seq= Fasta.read(inp)

    HHBLITS=HHBlits('/Users/tkaraman/opt/anaconda3/envs/python39/bin/hhblits', '/Users/tkaraman/hhsuite_dbs/UniRef30_2022_02_hhsuite/UniRef30_2022_02', n_iter=2, psi=True, hhm=True)
    HHBLITS.run(inp)
    hhm_mat=HHBLITS.process_hhm('out.hhm', fasta_string=seq)
    HHBLITS.plot_hhm(hhm_mat, outfile='HHM_map.pdf')
     
    CCMPRED=CCMPred('/Users/tkaraman/Programs/CCMpred/', n_top=500) 
    out=CCMPRED.run('./out.psi')
    CCMPRED.plotmat(out)
    CCMPRED.top_couplings(out, outfile='%s/%s_top_couplings_n500_s7.txt' %(path, prefix) ) 

    try:
        ids=Identifiers()
        ac       = ids.get_uniref_ac_from_hhr('%s/%s.hhr' %(path, prefix))
        pfam_id  = ids.get_pfam_id(ac)
        acs      = ids.get_all_uniref_acs_from_hhr('%s/%s.hhr' %(path, prefix))
        print(pfam_id)
        
        neff = ids.getNeff('%s/%s.hhr' %(path, prefix), cutoff=6)
        aligned_seqs=ids.get_aligned_seq_from_hhr('%s/%s.hhr' %(path, prefix), acs, seq)
       
    except: 
        print('NO HHR file for', prefix )
        continue

    #analysis
    nstructs=50
    outname='%s_couplings_analysis.txt'%prefix
    p = subprocess.run(["rm", "%s" %outname])

    for struct in list(aligned_seqs.keys())[:nstructs]:
        try:
            Pdb=PDB(struct)
            pdbfile      = Pdb.parse( struct, db='af' )
            coords       = Pdb.get_coords(pdbfile, atomname='[A-Z]+[1-9]{0,2}', chain='A', resid='[0-9]*')
            pdb_seq      = Pdb.get_seq(pdbfile, chain='A')
            offsets      = Pdb.get_offsets( pdb_seq, seq, pdb_resid_one=-1,  align=True)
    
            restr=restraints(struct)
            restr.calc_viol_rest('%s/%s_top_couplings_n500_s7.txt' %(path, prefix),  coords, n=20, dist_cutoff=5.1, conf_cutoff=0.0, offsets=offsets, outname=outname)           
            
        except RuntimeError as err:
            print(err.args)
            continue
     
    try:
        restr.analyse(outname)

    except: 
        continue


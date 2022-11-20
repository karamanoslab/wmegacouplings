import re
import numpy as np
import glob
import matplotlib.pyplot as plt
from pdb import PDB
from matplotlib.ticker import PercentFormatter

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['font.family'] = "sans-serif"


def analyse(analysis_file):

   #  pdbid     = analysis_file.split('/')[-1].split('_')[0]
#     Pdb       = PDB(pdbid)
#     homooligo = Pdb.is_homooligomer(pdbid)
    
    pattern  = re.compile('Fully violated restr ([0-9]*)')
    with open(analysis_file) as f:
        af = f.read()
    
    matches=re.findall(pattern, af)
    
    if len(matches) >12:
        raise RuntimeError('Too many violations in %s. Probably a short sequence; or offset calculation wrong' %analysis_file)
    
    dists,  scores,  pldtts=[], [], []
    vdists, vscores, vpldtts=[], [], []
    
    for viol in matches:
        patt=re.compile('# %s .*'%viol)
        viol_restrs=re.findall(patt, af)
                
        for viol_restr in viol_restrs:
            fields=viol_restr.split()
            coupling, pos1, pos2, plddt1, plddt2, dist, score =fields[1], fields[3], fields[5], float(fields[10]), float(fields[15]), fields[19], float(fields[21])
            mean_plddt=(float(plddt1) + float(plddt2))/2.
            
            if mean_plddt>95:# and score >1.0:
                print(viol_restr)
                vdists.append(float(dist))
                vscores.append(float(score))
                vpldtts.append( mean_plddt )
                
    
    vout=np.zeros( (len(vscores), 3))
    vout[:,0]=np.array(vdists)
    vout[:,1]=np.array(vpldtts)
    vout[:,2]=np.array(vscores)

    
    
    patt1=re.compile('# .*OK')
    all_restrs=re.findall(patt1, af)
    for restr in all_restrs:
        #print(restr)
        fields =restr.split()
        coupling, pos1, pos2, plddt1, plddt2, dist, score =fields[1], fields[3], fields[5], fields[10], fields[15], fields[19], fields[21]
        dists.append(float(dist))
        scores.append(float(score))
        pldtts.append( (float(plddt1) + float(plddt2))/2. )
       
    
    out=np.zeros( (len(scores), 3))
    out[:,0]=np.array(dists)
    out[:,1]=np.array(pldtts)
    out[:,2]=np.array(scores)

    
    f.close()
    return out, vout
          
# allfiles=[]
# for i in ['2', '3', '4', '5', '6', '7']:
#     path='./ccmpred%s_analysis_n50/'%i
#     files=glob.glob('%s/*_couplings_analysis.txt' %path)
#     allfiles += files 

allfiles=glob.glob('./ccmpred_all/*_couplings_analysis.txt')
_ok, _viol=np.empty( (1 ,3) ), np.empty( (1 ,3) )

for f in allfiles:
   print(f)
   try:
       a=analyse(f)
       _ok=np.concatenate( (_ok, a[0]), axis=0) 
       _viol=np.concatenate( (_viol, a[1]), axis=0) 
       
   except RuntimeError as err:
       print(err.args)


bins=100
fig=plt.figure(12, figsize=(10,3))
plt.rc('font', size =12)
ax=plt.subplot(131)
f=ax.hist(_ok[:,0], color='b', density=True, stacked=False, histtype='bar', alpha=0.5, range=(0,50),  bins=bins)
l=ax.hist(_viol[:,0],color='r', density=True, stacked=False, histtype='bar', alpha=0.5, range=(0,50),  bins=bins)
#ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.plot([5.1, 5.1], [0, 0.25], '--k')
ax.set_xlabel('Distance')
#ax.set_ylabel('#Occurences')
plt.xticks([0, 10, 20, 30, 40, 50])
print('Mean dist non viol: %.3f +/- %.3f' %(np.mean(_ok[:,0]), np.std(_ok[:,0]) ))
print('Mean dist viol: %.3f +/- %.3f\n' %(np.mean(_viol[:,0]), np.std(_viol[:,0]) ))

ax=plt.subplot(132)
ax.hist(_ok[:,1], color='b', density=True,  stacked=False,  alpha=0.5, range=(40,100), bins=bins)
ax.hist(_viol[:,1], color='r', density=True, stacked=True, alpha=0.5,  range=(40,100), bins=bins)
#ax.set_ylim([0.0,0.05])
ax.set_xlabel('mean pLDDT')
#ax.set_ylabel('#Occurences')
print('Mean pLDDT non viol: %.3f +/- %.3f' %(np.mean(_ok[:,1]), np.std(_ok[:,1]) ))
print('Mean pLDDT viol: %.3f +/- %.3f\n' %(np.mean(_viol[:,1]), np.std(_viol[:,1]) ))
        
ax=plt.subplot(133)
ax.hist(_ok[:,2], color='b', density=True,  stacked=True,  alpha=0.5, bins=bins)
ax.hist(_viol[:,2], color='r', density=True, stacked=True, alpha=0.5, bins=bins)
ax.set_xlabel('CCMpred conf. score')
#ax.set_ylabel('#Occurences')
print('Mean score non viol: %.3f +/- %.3f' %(np.mean(_ok[:,2]), np.std(_ok[:,2]) ))
print('Mean score viol: %.3f +/- %.3f' %(np.mean(_viol[:,2]), np.std(_viol[:,2]) ))

plt.tight_layout()
plt.savefig('distributions.pdf', transparent=True)
plt.show()

fig=plt.figure(13, figsize=(8,8))
plt.rc('font', size =12)
print(np.corrcoef(_viol[:,0], _viol[:,1])[0][1])
plt.plot(_viol[:,0], _viol[:,1], 'o', markersize=3, markerfacecolor='None', markeredgecolor='r' )
plt.show()



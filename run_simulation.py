import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import matplotlib.pyplot as plt
from data.simulate import Simulate
from tools.ccmpred import CCMPred



fasta= '''>test
KNPDDVFREFEKNPDDVFREFEKNPDDVFREFE'''
a='\
-----------|-----|--||---------------------------------'


coev1= {('W', 'Y'):(17, 20), ('F', 'Y'):(17, 20)}
coev2= {('K', 'E'):(11, 21), ('R', 'D'):(11, 21),  ('Q', 'N'):(11, 21)}  

print(fasta, file=open("sim.fasta", 'w') )


sim=Simulate()

nsteps=1
nmutants=1500

for step  in range(nsteps):
    np.random.seed()
    mutants, mutated_seqs=[], []
    outfile=open('sim.psi', 'w')
    outfile=open('sim.psi', 'a')
    for seq_record in SeqIO.parse("sim.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
    
        mutable_seq = MutableSeq(seq_record)

        for i in range(nmutants):
            mutated_seq=sim.mutate(mutable_seq, coev1, coev2, high_entropy_cols=[10, 30], low_entropy_cols=[],\
                                   coupling_strenghts=(10, 4, 4), coupled=True, mutation_rate=10.)
            print(i, mutated_seq)
            mutant=SeqRecord(Seq(mutated_seq), id="%i"%i)
            mutated_seqs.append(mutated_seq)
            print('Mutant%i %s' %(i,mutated_seq), file=outfile )
            
                               
    codes=sim.integer_encoding(mutated_seqs)  

    outfile.close()


    CCMPRED=CCMPred('/Users/tkaraman/Programs/CCMpred/') 
    out=CCMPRED.run('./sim.psi')
    plt.figure(10) 
    CCMPRED.plotmat(out)
    scores=CCMPRED.get_scores(out,[  (11, 20), 
                              (11, 21),
                              (17, 21),
                              (11, 20)  ])


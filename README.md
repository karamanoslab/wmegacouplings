  # WmegaCouplings

1. Install HHsuite (easy to do in anaconda)
2. Install CCMpred (cmake -DWITH_CUDA=OFF ${PATH_TO_CCMPRED_REPO to isntall without a gpu}
3. make sure that the CCMpred is in your PATH
4. Needs an internet connection
5. Written in python3
6. Dependencies numpy, matplotlib, biopython, seaborn
7. Run as 
   python run_for_af.py path_to_Uniref_db input.fasta
8. The x_coupling _analysis.txt contains the results :

   restr_number fasta_residue1(numbering starts from 0)  fasta_residue2  pdb_res1 atom1 plddt1 pdb_res2 atom2 plddt2_ accession_number distance, ccmpred_score
   
   At the end of the file you can find the pairs that are always violated
   
Output.zip contains the results of the analysis for 2012 Pfam families

Work in progress

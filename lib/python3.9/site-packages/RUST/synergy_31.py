
import numpy as np
from sys import argv
import commands, tempfile, shutil, os, sys

amino_acids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

tmp_dir = tempfile.mkdtemp( dir='.' )

tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_aligns_file_name = tmp_aligns_file.name
tmp_aligns_file.close()


infileopen = open("%s"%str(argv[2]) )
infileopen.readline()
list_amino = []
list_zscores = []
list_fold_change = []
list_loc  = []

for line in infileopen :
    linesplit = line[:-1].split(",")
    if len(linesplit) == 1 :break
    amino = linesplit[0]
    coverage =  map(float, linesplit[1:])
    coverage_a = coverage[0]
    if coverage_a == 0: continue
    coverage_n = [n/coverage_a for n in coverage[1:]]

    if len(amino) != 3 or len(set(amino) - set(amino_acids)) != 0 :
        sys.stderr.write( 'Tripeptide metafootprint file not in correct, check input file'  )
        if os.path.exists( tmp_dir ): shutil.rmtree( tmp_dir )
        exit()
    aminoA = amino[0]
    aminoB = amino[1]
    aminoC = amino[2]
        

    infileopen2 = open("%s"%str(argv[1]) )
    infileopen2.seek(0)
    infileopen2.readline()
    for line2 in infileopen2 :
        linesplit = line2[:-1].split(",")
        if len(linesplit) == 1 :break
        amino2 = linesplit[0]
        if len(amino2) != 1 or len(set(amino2) - set(amino_acids)) != 0 :
            sys.stderr.write( 'Amino acid metafootprint file not correct, check input file'  )
            if os.path.exists( tmp_dir ): shutil.rmtree( tmp_dir )
            exit()
        if amino2 in amino:
            coverage =  map(float, linesplit[1:])
            coverage_a = coverage[0]
            if coverage_a == 0: continue
            if amino2 == aminoA:
                coverage_n1 = [n/coverage_a for n in coverage[1:]]
            if amino2 == aminoB:
                coverage_n2 = [n/coverage_a for n in coverage[1:]]
            if amino2 == aminoC:
                coverage_n3 = [n/coverage_a for n in coverage[1:]]
    infileopen2.close()

    coverage_n_e = 0
    differences = []

    for number_i in range(11):
        coverage_n_e = coverage_n1[number_i]*coverage_n2[number_i+1]*coverage_n3[number_i+2]
        differences.append(abs(coverage_n[number_i]) - abs(coverage_n_e))
    for number_i in range(47,58):
        coverage_n_e = coverage_n1[number_i]*coverage_n2[number_i+1]*coverage_n3[number_i+2]
        differences.append(abs(coverage_n[number_i]) - abs(coverage_n_e))
        
    std_diff = np.std(differences)

    line_count = 0
    for number_i in range(0,len(coverage_n)-2):
        coverage_n_e = coverage_n1[number_i]*coverage_n2[number_i+1]*coverage_n3[number_i+2]
        
        list_amino.append(amino)
        list_zscores.append((coverage_n[number_i] - coverage_n_e)/std_diff)
        list_loc.append(number_i)
        if coverage_n_e == 0: 
            list_fold_change.append("not defined")
        else :
            list_fold_change.append(coverage_n[number_i] /coverage_n_e)
        
outfile = open("%s"%tmp_aligns_file_name,"w")
outfile.write("Tripeptide, Standard score, distance of 1st residue from A-site, fold change\n")
zipped_list = zip(list_zscores, list_amino, list_loc,list_fold_change)
zipped_list.sort()
zipped_list.reverse()
for zscore, amino, loc,fold_change in zipped_list :
    if abs(zscore) > 5 :
        outfile.write("%s, %s, %s, %s\n"%(amino,zscore, loc-40,fold_change))
            
outfile.close()		
shutil.move( tmp_aligns_file_name,str(argv[3]) )

        
if os.path.exists( tmp_dir ):
    shutil.rmtree( tmp_dir )
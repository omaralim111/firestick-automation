import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from pylab import log2, MaxNLocator
from sys import argv
from Bio import SeqIO
import numpy as np
import commands, tempfile, shutil, os, pysam, re, sys

def stop_err( msg, tmp_dir2):
    sys.stderr.write( '%s\n' % msg )
    if os.path.exists( tmp_dir2 ):
        shutil.rmtree( tmp_dir2 )
    sys.exit()
    
    
def RUST_metaprofile(infileopen36,ax36) :
    infileopen36.seek(0)
    infileopen36.readline()
    while 1:
        line = infileopen36.readline()
        linesplit = line.split(",")
        if len(linesplit) == 1 :break
        nucleotide_type = linesplit[0]
        coverage =  map(float, linesplit[1:])
        coverage_a = coverage[0]
        if coverage_a == 0: continue
        coverage_n = [n/coverage_a for n in coverage[1:]]
        ax36.plot(coverage_n[:-2],color = "gray")
      
    ax36.set_xticks([5,10,15,20,25,30,35,40,45,50,55])
    ax36.set_xticklabels([-35,-30,-25,-20,-15,-10,-5,0,5,10,15 ])
    ax36.set_xlabel("distance from A-site [codon]")
    ax36.set_ylabel("Tripeptide RUST ratio (observed/expected)")
    ax36.axvline(40, color = "red")
    #ax36.legend()
    

mpl.rcParams["xtick.direction"]        = "out"
mpl.rcParams["ytick.direction"]        = "out"
mpl.rcParams["legend.fontsize"]             = 10
mpl.rcParams["ytick.labelsize"]             = 10
mpl.rcParams["xtick.labelsize"]             = 10
mpl.rcParams['font.size']                   = 10
mpl.rcParams["axes.titlesize"]     = 10
mpl.rcParams["legend.frameon"]              = 0
mpl.rcParams["axes.axisbelow"]      = False
mpl.rcParams["xtick.major.pad"]     = 2.
mpl.rcParams["ytick.major.pad"]     = 2
mpl.rcParams["xtick.major.size"]    = 2.
mpl.rcParams["ytick.major.size"]     = 2
mpl.rcParams['axes.linewidth']     = 0.5
mpl.rcParams['ytick.major.width']     = 0.25
mpl.rcParams['xtick.major.width']     = 0.25
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams["legend.borderpad"]    = 0.01
mpl.rcParams["legend.labelspacing"] = 0.05
mpl.rcParams["legend.columnspacing"]= 0.5
mpl.rcParams["legend.borderaxespad"]=0.15
mpl.rcParams["legend.handlelength"]    = 1

tmp_dir = tempfile.mkdtemp( dir='.' )

tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_aligns_file_name = tmp_aligns_file.name
tmp_aligns_file.close()

sorted_bam_file = '%s.bam' % tmp_aligns_file_name
os.symlink(str(argv[2]), sorted_bam_file )
command = "samtools index %s"%(sorted_bam_file)
commands.getstatusoutput(command)

mRNA_sequences          =  str(argv[1])   #path to fastq file of transcripts
in_seq_handle           = open(mRNA_sequences)
cds_start_dict  = {}
cds_end_dict    = {}
for line in in_seq_handle:
    if line[0] != ">" : continue
    try :
        transcript_split= line.split(",")
        transcript       = transcript_split[0][1:]
        cds_start_dict[transcript]       = int(transcript_split[1])
        cds_end_dict[transcript]         = int(transcript_split[2])
    except: pass
in_seq_handle.seek(0)
seq_dict                = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
in_seq_handle.close()

offset          = int(argv[3])   #path to fastq file of transcripts
readlen_range   = str(argv[4])
readlen_rangesplit      = readlen_range.split(":")
if len(readlen_rangesplit) == 1: 
    accepted_read_lengths = [int(readlen_rangesplit[0])]
elif len(readlen_rangesplit) == 2 :
    accepted_read_lengths = [readlen for readlen in range(int(readlen_rangesplit[0]),int(readlen_rangesplit[1])+1)]
else : 
    stop_err( "Lengths of footprints parameter not in correct format, it should be either colon seperated with the second value greater or equal to the first, (28:32) or a single interger (31)",tmp_dir)
if len(accepted_read_lengths) == 0  :
    stop_err( "Lengths of footprints parameter not in correct format, it should be either colon seperated with the second value greater or equal to the first, (28:32) or a single interger (31)",tmp_dir)
    
amino_acids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
aligments_A1    = pysam.Samfile("%s"%sorted_bam_file, 'rb')  #path to aligments in bam format

tripeptide_enrichment_dict = {}
codon_enrichment_expected_dict= {}
for amino_acid in amino_acids :
    for amino_acid2 in amino_acids :
        for amino_acid3 in amino_acids :
            tripeptide   = "%s%s%s"%(amino_acid,amino_acid2,amino_acid3)
            tripeptide_enrichment_dict[tripeptide] = {}
            codon_enrichment_expected_dict[tripeptide] = []
            for number in range(0,60,1) :
                tripeptide_enrichment_dict[tripeptide][number] = [0.0,0.0]
                    


list_transcripts = seq_dict.keys()
for transcript in list_transcripts :
    try: #use supplied CDS annotation
        cds_start       = cds_start_dict[transcript]
        cds_end         = cds_end_dict[transcript]
        if cds_end      < cds_start: raise Exception
    except Exception:  #find longest ORF
        transcript_seq = str(seq_dict[transcript].seq)
        cds_start = -1
        start_post = []
        end_post = []        
        for match in re.finditer(r'(?=(%s))' % re.escape("ATG"), transcript_seq): start_post.append(match.start())
        for match in re.finditer(r'(?=(%s))' % re.escape("TAG" ),transcript_seq ):end_post.append(match.start())
        for match in re.finditer(r'(?=(%s))' % re.escape("TAA" ),transcript_seq ): end_post.append(match.start())
        for match in re.finditer(r'(?=(%s))' % re.escape("TGA" ),transcript_seq ):end_post.append(match.start())
        
        end_post.sort()
        len_max_orf = 0
        for value in start_post:
            for value2 in end_post:
                if value < value2:
                    if value%3 == value2%3:
                            len_orf = (value2-value)
                            if len_orf > len_max_orf:
                                cds_start = value
                                cds_end = value2+3
                                len_max_orf = len_orf
                            break
        if cds_start == -1:
            sys.stdout.write( '%s, AUG codon not found\n'%transcript  )
            continue
       
       
    elongation_region_all  = seq_dict[transcript][cds_start:cds_end].seq  
    elongation_region_part = elongation_region_all[120:-60]  # first 120 and last 60 nt are not used
    if len( elongation_region_part)% 3 != 0 : 
        sys.stdout.write( '%s, CDS not divisible by 3\n'%transcript  )
        continue
    peptide_sequence       = elongation_region_all.translate()

    profile_list            = [0.0 for n in range(cds_start+120, cds_end-60)]  # records ribo-seq profile
    if len(profile_list)    < 50: 
        sys.stdout.write( '%s, ORF too short\n'%transcript  )
        continue
    all_reads   = aligments_A1.fetch(transcript)

    len_elongation_region = len(profile_list)
    for read in all_reads :
        readlen = read.qlen
        if readlen not in accepted_read_lengths : continue          # selection of read of acceptable length
        A_site = read.pos+offset-cds_start-120                      # addition of offset
        if len_elongation_region > A_site > -1 :
            profile_list[A_site] += 1            
    average_gene_density = float(sum(profile_list))/len(profile_list)  #average gene density calculated

    if average_gene_density != 0 :
        num_codon = len([1 for number88 in range(0,len(profile_list),3) if ((profile_list[number88]+profile_list[number88+1]+profile_list[number88+2])/3)>average_gene_density])  
                                # number of codons that exceed average gene density 
        expected_codon_density = float(num_codon)/(len(profile_list)/3)    # expected enrichment value
        
        peptide_start 		= 0
        for sliding_w_n in range(0,len( elongation_region_part),3) :    #sliding window using increments of 3 nts
            amino_window = str(peptide_sequence[peptide_start:peptide_start+60])
            if len(set(amino_window) - set(amino_acids)) != 0: 
                peptide_start	+= 1
                continue

            for number in range(0,len(amino_window)-2) :
                amino_acid_3 = amino_window[number:number+3]
                tripeptide_enrichment_dict[amino_acid_3][number][0] += 1
                if (profile_list[sliding_w_n]+profile_list[sliding_w_n+1]+profile_list[sliding_w_n+2])/3 > average_gene_density:
                    tripeptide_enrichment_dict[amino_acid_3][number][1] += 1
            
            amino_acid_3 = amino_window[40:43]
            codon_enrichment_expected_dict[amino_acid_3].append(expected_codon_density)
            peptide_start += 1

tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_rust_file_name = tmp_aligns_file.name
tmp_aligns_file.close()

outfile = open("%s"%tmp_rust_file_name, "w")
outfile.write("tripeptide, expected value")
for number106 in range(-40,20): outfile.write(", %s"%number106)
outfile.write("\n")

list_codons = []
list_amino_acids 		= tripeptide_enrichment_dict.keys()
list_amino_acids.sort()
for amino2 in list_amino_acids :
    if amino2 in list_codons :continue
    list_codons.append(amino2)
    outfile.write("%s"%amino2)
    if codon_enrichment_expected_dict[amino2] != []:
        outfile.write( ", %s"%np.mean(codon_enrichment_expected_dict[amino2]))
    
    for number in range(0,60) :
        if tripeptide_enrichment_dict[amino2][number][0]  != 0 :
            outfile.write(  ", %s"%(tripeptide_enrichment_dict[amino2][number][1]/tripeptide_enrichment_dict[amino2][number][0]))
        else :
            outfile.write(  ", 0")
    outfile.write(  "\n")
outfile.close()

if not os.path.exists(str(argv[5]) ):
    os.mkdir(str(argv[5]) )
shutil.move( tmp_rust_file_name,str(argv[7]) )

try:
    fig = plt.figure(figsize = (6.69,6.0))  
    infileopen = open(str(argv[7]))
    ax1_metafootprint= fig.add_subplot(111)   
    RUST_metaprofile(infileopen, ax1_metafootprint)
    tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
    tmp_fig_file_name = tmp_aligns_file.name
    tmp_aligns_file.close()
    tmp_fig_file_name_png = "%s.png"%tmp_fig_file_name
    plt.savefig("%s"%tmp_fig_file_name_png)
    shutil.move( tmp_fig_file_name_png,os.path.join(str(argv[5]), 'RUST tripeptide metafootprint.png') )
    tmp_fig_file_name_png = "%s.svg"%tmp_fig_file_name
    plt.savefig("%s"%tmp_fig_file_name_png)
    shutil.move( tmp_fig_file_name_png,os.path.join(str(argv[5]), 'RUST tripeptide metafootprint.svg') )
    plt.clf()

        
    html = '<p><img alt="RUST tripeptide metafootprint" src="RUST tripeptide metafootprint.png" border="1"><br></p>'

    with open(str(argv[6]), 'w') as g:
        g.write(html)    
except:
    sys.stdout.write("Error producing images")
        
if os.path.exists( tmp_dir ):
    shutil.rmtree( tmp_dir )
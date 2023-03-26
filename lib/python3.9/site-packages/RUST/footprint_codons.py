from Bio import SeqIO
from sys import argv
import tempfile, shutil, os, commands, pysam,sys, re, numpy

def rank(lsit1):
    lsit2 = []
    lsit1s = lsit1[:]
    lsit1s.sort()
    dict_ranks = {}
    for value_i, value in enumerate(lsit1s):
        dict_ranks.setdefault(value,[]).append(value_i)
    
    for value in lsit1:
        lsit2.append(numpy.mean(dict_ranks[value]))
    return lsit2

def stop_err( msg, tmp_dir2):
    sys.stderr.write( '%s\n' % msg )
    if os.path.exists( tmp_dir2 ):
        shutil.rmtree( tmp_dir2 )
    sys.exit()
    
    

tmp_dir = tempfile.mkdtemp( dir='.' )
tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_aligns_file_name = tmp_aligns_file.name
tmp_aligns_file.close()
    
RUST_file = open(str(argv[5]))  # file output of RUST_script.py
RUST_file.readline()
codon_rust_dict = {}
for line in RUST_file :
    linesplit = line.split(",")
    if len(linesplit) == 1 :break
    codon = linesplit[0]
    if len(codon) != 3 or len(set(codon)-set(["A","T","G","C"]))!= 0:
        stop_err("Codon metafootprint file not correct, check input file",tmp_dir)
    codon_rust_dict[codon] = {}
    rust_values =  map(float, linesplit[1:])
    expected = rust_values[0]
    rust_metafootprint = [ro_value/expected for ro_value in rust_values[1:]]
    for n in range(34,46) :
        codon_rust_dict[codon][n-40] = rust_metafootprint[n]    #for 12 codons positions near A-site 
RUST_file.close()


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
        transcript_split= line.split("\t")
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
aligments_A1    = pysam.Samfile("%s"%sorted_bam_file, 'rb')


tmp_aligns_file3 = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_aligns_file3_name = tmp_aligns_file3.name
tmp_aligns_file3.close()
correlations_file = open("%s"%tmp_aligns_file3_name, "w")
correlations_file.write("transcript,average read density,Spearman's coefficient,Pearson's coefficient\n")


list_transcripts = seq_dict.keys()
for transcript in list_transcripts :
    try:
        cds_start       = cds_start_dict[transcript]
        cds_end         = cds_end_dict[transcript]
        if cds_end      < cds_start: raise Exception
    except Exception:
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
        
    elongation_region_all 	= str(seq_dict[transcript][cds_start:cds_end].seq )
    
    if len(elongation_region_all)% 3 != 0 :  	#genes with codon region not divisible by 3 skipped
        sys.stdout.write( '%s, CDS not divisible by 3\n'%transcript  )
        continue
    peptide_sequence        = seq_dict[transcript][cds_start:cds_end].seq.translate()

    profile_expect 	= []
    for n in range(0,len(elongation_region_all[120:-60]),3) :  # predicts profile from 120 nts after start to 60 before stop
        minus6_plus5_footprint = elongation_region_all[120+n-18:120+n+19] #contains sequence of region used to predict profile
        value = 1.0
        amino_loc = -6
        for number in range(0,len(minus6_plus5_footprint)-2,3) :
            codon = minus6_plus5_footprint[number:number+3]
            if len(set(codon) - set(["A","T","G","C"])) != 0 or codon in ["TAG","TGA","TAA"]: 
                amino_loc+= 1
                continue
            value = value* codon_rust_dict[codon][amino_loc]
            amino_loc+= 1
        profile_expect.append(value)
    profile_expect_sum = sum(profile_expect)
    profile_expect_probablility = [float(value)/profile_expect_sum for value in profile_expect]
    

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
    if average_gene_density > 0 :
        profiles_control_codon = [profile_list[codon_ind]+profile_list[codon_ind+1]+profile_list[codon_ind+2] for codon_ind in range(0,len(profile_list),3)]
        spearmanr_value = numpy.corrcoef(rank(profiles_control_codon), rank(profile_expect))[0, 1]
        pearsonr_value = numpy.corrcoef(profiles_control_codon, profile_expect)[0, 1]
        correlations_file.write("%s,%s,%s,%s\n"%(transcript,average_gene_density,spearmanr_value,pearsonr_value))
  

correlations_file.close()
shutil.move( tmp_aligns_file3_name,str(argv[6]) )

if os.path.exists( tmp_dir ):
    shutil.rmtree( tmp_dir )
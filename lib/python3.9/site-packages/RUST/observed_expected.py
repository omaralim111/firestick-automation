import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from sys import argv
import tempfile, shutil, os, commands, pysam, sys, numpy, re
from pylab import MaxNLocator

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
        codon_rust_dict[codon][n-40] = rust_metafootprint[n]   #for 12 codons positions near A-site the RUST ratios are recorded
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

tmp_aligns_file2 = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_aligns_file2_name = tmp_aligns_file2.name
tmp_aligns_file2.close()
open_file = open("%s"%tmp_aligns_file2_name, "w")


list_transcripts = seq_dict.keys()
transcript_of_inter = str(argv[6])
transcript_of_inter2 = transcript_of_inter
if transcript_of_inter not in list_transcripts:
    count_occurences = 0 
    for known_transcript in list_transcripts:
        if transcript_of_inter in known_transcript:
            transcript_of_inter2 = known_transcript
            count_occurences += 1
    if transcript_of_inter2 == transcript_of_inter:
        stop_err( "Transcript not in Transcriptome file",tmp_dir)
    if count_occurences > 1 :
        stop_err( "%s not unique identifier"%transcript_of_inter,tmp_dir)
    sys.stdout.write( '%s not in Transcriptome file, data provided for %s' %(transcript_of_inter,transcript_of_inter2 ))

for transcript in [transcript_of_inter2] :
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
        if cds_start == -1: continue
        
    elongation_region_all 	= str(seq_dict[transcript][cds_start:cds_end].seq )
    if len(elongation_region_all)% 3 != 0 :     #genes with codon region not divisible by 3 skipped
        stop_err( '%s, CDS not divisible by 3\n'%transcript ,tmp_dir )

    profile_expect 	= []
    for n in range(0,len(elongation_region_all[120:-60]),3) :  # predicts profile from 120 nts after start to 60 before stop
        minus6_plus5_footprint = elongation_region_all[120+n-18:120+n+16] #contains sequence of region used to predict profile
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
    all_reads   = aligments_A1.fetch(transcript)

    len_elongation_region = len(profile_list)
    for read in all_reads :
        readlen = read.qlen
        if readlen not in accepted_read_lengths : continue          # selection of read of acceptable length
        A_site = read.pos+offset-cds_start-120                      # addition of offset
        if len_elongation_region > A_site > -1 :
            profile_list[A_site] += 1
            
    profiles_control_codon = [profile_list[codon_ind]+profile_list[codon_ind+1]+profile_list[codon_ind+2] for codon_ind in range(0,len(profile_list),3)]
    profile_expect_probablility_index = 0
    open_file.write("%s\n"%transcript)
    open_file.write("codon, expected probability, alignments\n")
    for coordinate_index in range(0,len(elongation_region_all[120:-60]),3) :  
        codon = elongation_region_all[120+coordinate_index:120+coordinate_index+3]
        open_file.write("%s, "%(codon))
        open_file.write("%s, "%(profile_expect_probablility[profile_expect_probablility_index] ))
        open_file.write("%s\n"%(profiles_control_codon[profile_expect_probablility_index] ))
        profile_expect_probablility_index += 1
open_file.close()
if not os.path.exists(str(argv[7]) ):
    os.mkdir(str(argv[7])) 

str_i =     'profile_data_%s.csv'%transcript
shutil.move( tmp_aligns_file2_name,os.path.join(str(argv[7]),str_i ) )

    
fig = plt.figure(figsize = (6.69,6.0))  
plt.subplots_adjust(left =0.09, right=.87)
ax = fig.add_subplot(111)
ax.plot(profiles_control_codon,color = "gray",label = "observed")
ax2 = ax.twinx()
ax2.plot(profile_expect_probablility,"--", color = "DarkMagenta",label = "predicted")

ax.text(0.1,1.05, "r =%s"%round(numpy.corrcoef(profiles_control_codon, profile_expect_probablility)[0, 1],2), transform=ax.transAxes)
    

l= ax.legend( bbox_to_anchor=(0, 0, 0.890, 1.05), bbox_transform=ax.transAxes, ncol =1)
l= ax2.legend( bbox_to_anchor=(0, 0, 0.890, 1.10), bbox_transform=ax2.transAxes, ncol =1)

ax.set_xlabel("transcript coordinates [codon]")
ax.set_ylabel("# alignments")
ax.yaxis.set_major_locator(MaxNLocator(5))
tciks = ax.get_xticks()
cds_start_codon = cds_start/3
tciks2 = [int(n)+40+cds_start_codon for n in tciks]
ax.set_xticklabels(tciks2)

ax.set_title(transcript_of_inter)
for tl in ax2.get_yticklabels():
    tl.set_color('DarkMagenta')
ax2.yaxis.set_major_locator(MaxNLocator(5))
ax2.set_ylabel("probability", color = "darkmagenta")    
ax2.set_xlim(0, len(profile_expect_probablility))
ax.set_xlim(0, len(profile_expect_probablility))

tmp_aligns_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
tmp_fig_file_name = tmp_aligns_file.name
tmp_aligns_file.close()

tmp_fig_file_name_png = "%s.png"%tmp_fig_file_name
plt.savefig("%s"%tmp_fig_file_name_png)
shutil.move( tmp_fig_file_name_png,os.path.join(str(argv[7]), 'profile_%s.png'%transcript) )
tmp_fig_file_name_png = "%s.svg"%tmp_fig_file_name
plt.savefig("%s"%tmp_fig_file_name_png)
shutil.move( tmp_fig_file_name_png,os.path.join(str(argv[7]), 'profile_%s.svg'%transcript) )


html = '<p><img alt="Observed_expected_profile" src="profile_%s.png" border="1"><br><a href="profile_data_%s.csv">profile_data_%s.csv</a></p>'%(transcript,transcript,transcript)

with open(str(argv[8]), 'w') as g:
    g.write(html)    
    
if os.path.exists( tmp_dir ):
    shutil.rmtree( tmp_dir )
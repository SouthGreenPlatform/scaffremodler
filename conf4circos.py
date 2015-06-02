import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR):
	print cmd_line
	try:
		tmp = (tempfile.NamedTemporaryFile().name)+'.error'
		# print tmp
		error = open(tmp, 'w')
		proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
		returncode = proc.wait()
		error.close()
		error = open( tmp, 'rb' )
		stderr = ''
		buffsize = 1048576
		try:
			while True:
				stderr += error.read( buffsize )
				if not stderr or len( stderr ) % buffsize != 0:
					break
		except OverflowError:
			pass
		error.close()
		os.remove(tmp)
		if returncode != 0:
			raise Exception, stderr
	except Exception, e:
		stop_err( ERROR + str( e ) )

def create_kar(REF, CHR, OUT_KAR, OUT_N):
	record_dict = SeqIO.index(REF, "fasta")
	file = open(CHR)
	out1 = open(OUT_KAR,'w')
	for line in file:
		data = line.split()
		if data:
			out1.write("chr - "+data[0]+ " "+data[0]+" 0 "+str(len(str(record_dict[data[0]].seq))+1)+" blue\n")
	file.close()
	
	out2 = open(OUT_N,'w')
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			sequence = str(record_dict[data[0]].seq)
			pos = 0
			debut = 0
			sur_N = 0
			if sequence[0] == 'N' or sequence[0] == 'n':
				sys.exit('the program cannot work because sequence '+data[0]+' begins with N')
			elif sequence[-1] == 'N' or sequence[-1] == 'n':
				sys.exit('the program cannot work because sequence '+data[0]+' ends with N')
			out2.write(data[0]+'\t'+str(pos+1)+'\t'+str(pos+1)+'\t'+data[0]+'-'+str(pos+1)+'-'+'p\n')
			for n in sequence:
				if n == 'N' or n == 'n':
					if not(sur_N):
						sur_N = 1
						debut = pos
						out2.write(data[0]+'\t'+str(pos)+'\t'+str(pos)+'\t'+data[0]+'-'+str(pos)+'-'+'v\n')
				elif sur_N:
					sur_N = 0
					out1.write('band '+data[0]+' p36.33 p36.33 '+str(debut)+' '+str(pos+1)+' black\n')
					out2.write(data[0]+'\t'+str(pos+1)+'\t'+str(pos+1)+'\t'+data[0]+'-'+str(pos+1)+'-'+'p\n')
				pos += 1
			out2.write(data[0]+'\t'+str(pos)+'\t'+str(pos)+'\t'+data[0]+'-'+str(pos)+'-'+'f\n')
	out1.close()
	out2.close()
		
##########################################################################################################################################################
#fonction that create link file based on discordant zone
def create_discord_link(FILE, OUT):
	outfile = open(OUT,'w')
	i = 1
	file = open(FILE)
	for line in file:
		data = line.split()
		if data:
			if data[0][0] != "#":
				if data[13] == 'PASSED':
					if i == 1:
						outfile.write('link'+str(i)+' '+data[0]+' '+data[1]+' '+data[2]+'\n')
						outfile.write('link'+str(i)+' '+data[5]+' '+data[6]+' '+data[7])
						i = i + 1
					else:
						outfile.write('\nlink'+str(i)+' '+data[0]+' '+data[1]+' '+data[2])
						outfile.write('\nlink'+str(i)+' '+data[5]+' '+data[6]+' '+data[7])
						i = i + 1
	outfile.close()

##########################################################################################################################################################
#fonction that create paired-read link file based on .list file of ApMap
def create_read_link(FILE, OUT, TYPE):
	outfile = open(OUT,'w')
	i = 1
	file = open(FILE)
	for line in file:
		data = line.split()
		if data:
			if data[5] == 'discard':
				if data[6] == TYPE:
					if i == 1:
						outfile.write('link'+str(i)+' '+data[3]+' '+data[1]+' '+data[1]+'\n')
						outfile.write('link'+str(i)+' '+data[4]+' '+data[2]+' '+data[2])
						i = i + 1
					else:
						outfile.write('\nlink'+str(i)+' '+data[3]+' '+data[1]+' '+data[1])
						outfile.write('\nlink'+str(i)+' '+data[4]+' '+data[2]+' '+data[2])
						i = i + 1
			else:
				if data[5] == TYPE:
					if i == 1:
						outfile.write('link'+str(i)+' '+data[3]+' '+data[1]+' '+data[1]+'\n')
						outfile.write('link'+str(i)+' '+data[4]+' '+data[2]+' '+data[2])
						i = i + 1
					else:
						outfile.write('\nlink'+str(i)+' '+data[3]+' '+data[1]+' '+data[1])
						outfile.write('\nlink'+str(i)+' '+data[4]+' '+data[2]+' '+data[2])
						i = i + 1
	outfile.close()


##########################################################################################################################################################
#fonction that calculate coverage over a certain windows size on the reference

def mediane(LISTE):
	L = sorted(LISTE)
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return L[p]

def moyenne(L):
	if len(L) == 0:
		moyenne = 0
	else:
		moyenne = sum(L)/float(len(L))
	return moyenne


def calcul_couv_moy(COUV, CHR, WIN, OUT):
	#############################################
	#Recording chromosome size
	DIC = {}
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			DIC[data[0]] = int(data[1])
	file.close()
	#Done
	#############################################
	#calculating median of mean coverage
	LISTE_total = []
	FILE = open(COUV)
	chr = ""
	i = 0
	for LINE in FILE:
		DATA = LINE.split()
		if DATA != []:
			i = i + 1
			if chr == "":#on est au tout debut
				chr = DATA[0]
				LISTE = []
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						LISTE_total.append(sum(LISTE)/float(len(LISTE)))
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					LISTE_total.append(sum(LISTE)/float(len(LISTE)))
					LISTE = []
			elif chr != DATA[0]:#on arrive sur un nouveau chromosome
				while i <= DIC[chr]:
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						LISTE_total.append(sum(LISTE)/float(len(LISTE)))
						LISTE = []
					i = i + 1
				if len(LISTE)%WIN != 0:
					LISTE_total.append(sum(LISTE)/float(len(LISTE)))
				LISTE = []
				i = 1
				chr = DATA[0]
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						LISTE_total.append(sum(LISTE)/float(len(LISTE)))
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					LISTE_total.append(sum(LISTE)/float(len(LISTE)))
					LISTE = []
			else:
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						LISTE_total.append(sum(LISTE)/float(len(LISTE)))
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					LISTE_total.append(sum(LISTE)/float(len(LISTE)))
					LISTE = []
	MEDIANE = mediane(LISTE_total)
	#Done
	#############################################
	#calculating mean covearge
	LISTE_total = []
	OUTFILE = open(OUT,'w')
	FILE = open(COUV)
	chr = ""
	i = 0
	for LINE in FILE:
		DATA = LINE.split()
		if DATA != []:
			i = i + 1
			if chr == "":#on est au tout debut
				chr = DATA[0]
				LISTE = []
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
						LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
					LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
					LISTE = []
			elif chr != DATA[0]:#on arrive sur un nouveau chromosome
				while i <= DIC[chr]:
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
						LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
						LISTE = []
					i = i + 1
				if len(LISTE)%WIN != 0:
					OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
					LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
				LISTE = []
				i = 1
				chr = DATA[0]
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
						LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
					LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
					LISTE = []
			else:
				while i < int(DATA[1]):
					LISTE.append(0)
					if len(LISTE)%WIN == 0:
						OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
						LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
						LISTE = []
					i = i + 1
				LISTE.append(int(DATA[2]))
				if len(LISTE)%WIN == 0:
					OUTFILE.write(chr+'\t'+str(i-len(LISTE))+'\t'+str(i-1)+'\t'+str((sum(LISTE)/float(len(LISTE)))-MEDIANE)+'\n')
					LISTE_total.append((sum(LISTE)/float(len(LISTE)))-MEDIANE)
					LISTE = []
	OUTFILE.flush()
	OUTFILE.close()
	return [MEDIANE, moyenne(LISTE_total)]

def create_tile(FILE, OUT):
	outfile = open(OUT, 'w')
	file = open(FILE)
	mot = 'color=black'
	for line in file:
		data = line.split()
		if data:
			if data[0][0] != '#':
				if data[4] == 'W':
					outfile.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+mot+'\n')
					if mot == 'color=black':
						mot = 'color=yellow'
					else:
						mot = 'color=black'
	file.close()
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This program takes in input all files needed to generate circos and output severals files that will be used to generate different circos plus a config file")
	
	# Wrapper options.
	#For input
	
	parser.add_option( '', '--ref', dest='ref', default='not_filled', help='The multi-fasta reference file')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='File containing col1: chromosome name and col2: chromsome size')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	
	parser.add_option( '', '--cov', dest='cov', default='not_filled', help='The coverage file')
	parser.add_option( '', '--window', dest='window', default=1000, help='Window size (base pair), [default: %default]')
	
	parser.add_option( '', '--frf', dest='frf', default='not_filled', help='frf_discord.score file')
	parser.add_option( '', '--ff', dest='ff', default='not_filled', help='ff_discord.score file')
	parser.add_option( '', '--rr', dest='rr', default='not_filled', help='rr_discord.score file')
	parser.add_option( '', '--ins', dest='ins', default='not_filled', help='ins_discord.score file')
	parser.add_option( '', '--delet', dest='delet', default='not_filled', help='delet_discord.score file')
	parser.add_option( '', '--chr_rr', dest='chr_rr', default='not_filled', help='chr_rr_discord.score file')
	parser.add_option( '', '--chr_rf', dest='chr_rf', default='not_filled', help='chr_rf_discord.score file')
	parser.add_option( '', '--chr_ff', dest='chr_ff', default='not_filled', help='chr_ff_discord.score file')
	parser.add_option( '', '--chr_fr', dest='chr_fr', default='not_filled', help='chr_fr_discord.score file')
	
	parser.add_option( '', '--liste_read', dest='liste_read', default='not_filled', help='A file containing information on mapped reads (.list file of ApMap)')
	
	parser.add_option( '', '--dis_prop', dest='dis_prop', default='not_filled', help='The file containing proportions of discordant reads (.prop file of ApMap)')
	
	parser.add_option( '', '--agp', dest='agp', default='not_filled', help='An agp file locating scaffold along chromosomes')
	
	# For output
	parser.add_option( '', '--prefix', dest='prefix', default='not_filled', help='Prefix for all output files. If this options is passed, all others output options are ignored, [default: %default]')
	
	# If no prefix arguments are passed
	parser.add_option( '', '--out_kar', dest='out_kar', default='circos_karyotype.txt', help='Karyotype output file, [default: %default]')
	parser.add_option( '', '--out_N', dest='out_N', default='circos_loc_N.txt', help='File name of text file locating N, [default: %default]')
	
	parser.add_option( '', '--out_cov', dest='out_cov', default='circos_mean.cov', help='Mean coverage output file, [default: %default]')
	
	parser.add_option( '', '--out_frf', dest='out_frf', default='circos_zone_frf.link', help='A link output file name corresponding to frf_discord.score, [default: %default]')
	parser.add_option( '', '--out_ff', dest='out_ff', default='circos_zone_ff.link', help='A link output file name corresponding to ff_discord.score, [default: %default]')
	parser.add_option( '', '--out_rr', dest='out_rr', default='circos_zone_rr.link', help='A link output file name corresponding to rr_discord.score, [default: %default]')
	parser.add_option( '', '--out_ins', dest='out_ins', default='circos_zone_ins.link', help='A link output file name corresponding to ins_discord.score, [default: %default]')
	parser.add_option( '', '--out_delet', dest='out_delet', default='circos_zone_delet.link', help='A link output file name corresponding to delet_discord.score, [default: %default]')
	parser.add_option( '', '--out_chr_rr', dest='out_chr_rr', default='circos_zone_chr_rr.link', help='A link output file name corresponding to chr_rr_discord.score, [default: %default]')
	parser.add_option( '', '--out_chr_rf', dest='out_chr_rf', default='circos_zone_chr_rf.link', help='A link output file name corresponding to chr_rf_discord.score, [default: %default]')
	parser.add_option( '', '--out_chr_ff', dest='out_chr_ff', default='circos_zone_chr_ff.link', help='A link output file name corresponding to chr_ff_discord.score, [default: %default]')
	parser.add_option( '', '--out_chr_fr', dest='out_chr_fr', default='circos_zone_chr_fr.link', help='A link output file name corresponding to chr_fr_discord.score, [default: %default]')
		
	parser.add_option( '', '--Rout_rf', dest='Rout_rf', default='circos_read_rf.link', help='A link output file name corresponding to mapped fr reads, [default: %default]')
	parser.add_option( '', '--Rout_fr', dest='Rout_fr', default='circos_read_fr.link', help='A link output file name corresponding to mapped rf reads, [default: %default]')
	parser.add_option( '', '--Rout_ff', dest='Rout_ff', default='circos_read_ff.link', help='A link output file name corresponding to mapped ff reads, [default: %default]')
	parser.add_option( '', '--Rout_rr', dest='Rout_rr', default='circos_read_rr.link', help='A link output file name corresponding to mapped rr reads, [default: %default]')
	parser.add_option( '', '--Rout_ins', dest='Rout_ins', default='circos_read_ins.link', help='A link output file name corresponding to mapped ins reads, [default: %default]')
	parser.add_option( '', '--Rout_delet', dest='Rout_delet', default='circos_read_delet.link', help='A link output file name corresponding to mapped delet reads, [default: %default]')
	parser.add_option( '', '--Rout_chr_rr', dest='Rout_chr_rr', default='circos_read_chr_rr.link', help='A link output file name corresponding to mapped chr_rr reads, [default: %default]')
	parser.add_option( '', '--Rout_chr_rf', dest='Rout_chr_rf', default='circos_read_chr_rf.link', help='A link output file name corresponding to mapped chr_rf reads, [default: %default]')
	parser.add_option( '', '--Rout_chr_ff', dest='Rout_chr_ff', default='circos_read_chr_ff.link', help='A link output file name corresponding to mapped chr_ff reads, [default: %default]')
	parser.add_option( '', '--Rout_chr_fr', dest='Rout_chr_fr', default='circos_read_chr_fr.link', help='A link output file name corresponding to mapped chr_fr reads, [default: %default]')
	
	parser.add_option( '', '--out_scaff', dest='out_scaff', default='circos_scaffold.tile', help='A tile output file name corresponding to scaffolds, [default: %default]')
	
	parser.add_option( '', '--output', dest='output', default='config_circos.conf', help='The output of the conf file, [default: %default]')
	
	
	
	(options, args) = parser.parse_args()
	
	if options.ref == 'not_filled':
		sys.exit('--ref argument is missing')
	if options.chr == 'not_filled':
		sys.exit('--chr argument is missing')
		
	if options.prefix != 'not_filled':
		create_kar(options.ref, options.chr, options.prefix+'_karyotype.txt', options.prefix+'_loc_N.txt')
	else:
		create_kar(options.ref, options.chr, options.out_kar, options.out_N)
	
	if options.cov != 'not_filled':
		if options.cov == 'not_filled':
			sys.exit('--cov argument is missing')
		if options.window == 'not_filled':
			sys.exit('--window argument is missing')
		if options.out_cov == 'not_filled':
			sys.exit('--out_cov argument is missing')
		if options.prefix != 'not_filled':
			Value = calcul_couv_moy(options.cov, options.chr, float(options.window), options.prefix+'_mean.cov')
		else:
			Value = calcul_couv_moy(options.cov, options.chr, float(options.window), options.out_cov)
	
	#For discordant link
	if options.frf != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.frf, options.prefix+'_zone_frf.link')
		else:
			create_discord_link(options.frf, options.out_frf)
			
	if options.ff != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.ff, options.prefix+'_zone_ff.link')
		else:
			create_discord_link(options.ff, options.out_ff)
	if options.rr != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.rr, options.prefix+'_zone_rr.link')
		else:
			create_discord_link(options.rr, options.out_rr)
	if options.ins != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.ins, options.prefix+'_zone_ins.link')
		else:
			create_discord_link(options.ins, options.out_ins)
	if options.delet != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.delet, options.prefix+'_zone_delet.link')
		else:
			create_discord_link(options.delet, options.out_delet)
	if options.chr_rr != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.chr_rr, options.prefix+'_zone_chr_rr.link')
		else:
			create_discord_link(options.chr_rr, options.out_chr_rr)
	if options.chr_rf != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.chr_rf, options.prefix+'_zone_chr_rf.link')
		else:
			create_discord_link(options.chr_rf, options.out_chr_rf)
	if options.chr_ff != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.chr_ff, options.prefix+'_zone_chr_ff.link')
		else:
			create_discord_link(options.chr_ff, options.out_chr_ff)
	if options.chr_fr != 'not_filled':
		if options.prefix != 'not_filled':
			create_discord_link(options.chr_fr, options.prefix+'_zone_chr_fr.link')
		else:
			create_discord_link(options.chr_fr, options.out_chr_fr)
	
	#For read link
	if options.liste_read != 'not_filled':
		if options.orient == 'rf':
			if options.prefix != 'not_filled':
				create_read_link(options.liste_read, options.prefix+'_read_rf.link', 'ok')
				create_read_link(options.liste_read, options.prefix+'_read_fr.link', 'fr')
			else:
				create_read_link(options.liste_read, options.Rout_rf, 'ok')
				create_read_link(options.liste_read, options.Rout_fr, 'fr')
		elif options.orient == 'fr':
			if options.prefix != 'not_filled':
				create_read_link(options.liste_read, options.prefix+'_read_rf.link', 'rf')
				create_read_link(options.liste_read, options.prefix+'_read_fr.link', 'ok')
			else:
				create_read_link(options.liste_read, options.Rout_rf, 'rf')
				create_read_link(options.liste_read, options.Rout_fr, 'ok')
		else:
			mot = 'Unrecognized argument in --orient %s' % options.orient
			sys.exit(mot)
		if options.prefix != 'not_filled':
			create_read_link(options.liste_read, options.prefix+'_read_ff.link', 'ff')
			create_read_link(options.liste_read, options.prefix+'_read_rr.link', 'rr')
			create_read_link(options.liste_read, options.prefix+'_read_ins.link', 'ins')
			create_read_link(options.liste_read, options.prefix+'_read_delet.link', 'del')
			create_read_link(options.liste_read, options.prefix+'_read_chr_rr.link', 'chr_rr')
			create_read_link(options.liste_read, options.prefix+'_read_chr_rf.link', 'chr_rf')
			create_read_link(options.liste_read, options.prefix+'_read_chr_ff.link', 'chr_ff')
			create_read_link(options.liste_read, options.prefix+'_read_chr_fr.link', 'chr_fr')
		else:
			create_read_link(options.liste_read, options.Rout_ff, 'ff')
			create_read_link(options.liste_read, options.Rout_rr, 'rr')
			create_read_link(options.liste_read, options.Rout_ins, 'ins')
			create_read_link(options.liste_read, options.Rout_delet, 'del')
			create_read_link(options.liste_read, options.Rout_chr_rr, 'chr_rr')
			create_read_link(options.liste_read, options.Rout_chr_rf, 'chr_rf')
			create_read_link(options.liste_read, options.Rout_chr_ff, 'chr_ff')
			create_read_link(options.liste_read, options.Rout_chr_fr, 'chr_fr')
		
	if options.agp != 'not_filled':
		if options.prefix != 'not_filled':
			create_tile(options.agp, options.prefix+'_scaffold.tile')
		else:
			create_tile(options.agp, options.out_scaff)

	config = ConfigParser.RawConfigParser()
	if options.prefix != 'not_filled':
		config.add_section('General')
		config.set('General','chr', options.chr)
		config.set('General','out_kar', options.prefix+'_karyotype.txt')
		config.set('General','out_N', options.prefix+'_loc_N.txt')
		config.set('General','orient', options.orient)
		if options.cov != 'not_filled':
			config.add_section('Coverage')
			config.set('Coverage','cov', options.prefix+'_mean.cov')
			config.set('General','cov', 'yes')
			config.set('Coverage','median_cov', str(Value[0]))
			config.set('Coverage','mean_cov', str(Value[1]))
		else:
			config.set('General','cov', 'no')
		config.add_section('Discord_link')
		config.add_section('Discord_zone')
		if options.frf != 'not_filled':
			config.set('Discord_link','frf', options.prefix+'_zone_frf.link')
			config.set('General','frf', 'yes')
		else:
			config.set('General','frf', 'no')
			
		if options.ff != 'not_filled':
			config.set('Discord_link','ff', options.prefix+'_zone_ff.link')
			config.set('General','ff', 'yes')
		else:
			config.set('General','ff', 'no')
			
		if options.rr != 'not_filled':
			config.set('Discord_link','rr', options.prefix+'_zone_rr.link')
			config.set('General','rr', 'yes')
		else:
			config.set('General','rr', 'no')
			
		if options.ins != 'not_filled':
			config.set('Discord_link','ins', options.prefix+'_zone_ins.link')
			config.set('General','ins', 'yes')
		else:
			config.set('General','ins', 'no')
			
		if options.delet != 'not_filled':
			config.set('Discord_link','delet', options.prefix+'_zone_delet.link')
			config.set('General','delet', 'yes')
		else:
			config.set('General','delet', 'no')
			
		if options.chr_rr != 'not_filled':
			config.set('Discord_link','chr_rr', options.prefix+'_zone_chr_rr.link')
			config.set('Discord_zone','chr_rr', options.chr_rr)
			config.set('General','chr_rr', 'yes')
		else:
			config.set('General','chr_rr', 'no')
			
		if options.chr_rf != 'not_filled':
			config.set('Discord_link','chr_rf', options.prefix+'_zone_chr_rf.link')
			config.set('Discord_zone','chr_rf', options.chr_rf)
			config.set('General','chr_rf', 'yes')
		else:
			config.set('General','chr_rf', 'no')
			
		if options.chr_fr != 'not_filled':
			config.set('Discord_link','chr_fr', options.prefix+'_zone_chr_fr.link')
			config.set('Discord_zone','chr_fr', options.chr_fr)
			config.set('General','chr_fr', 'yes')
		else:
			config.set('General','chr_fr', 'no')
			
		if options.chr_ff != 'not_filled':
			config.set('Discord_link','chr_ff', options.prefix+'_zone_chr_ff.link')
			config.set('Discord_zone','chr_ff', options.chr_ff)
			config.set('General','chr_ff', 'yes')
		else:
			config.set('General','chr_ff', 'no')
			
		if options.liste_read != 'not_filled':
			config.add_section('Read_link')
			config.set('Read_link','rf', options.prefix+'_read_rf.link')
			config.set('Read_link','fr', options.prefix+'_read_fr.link')
			config.set('Read_link','ff', options.prefix+'_read_ff.link')
			config.set('Read_link','rr', options.prefix+'_read_rr.link')
			config.set('Read_link','ins', options.prefix+'_read_ins.link')
			config.set('Read_link','del', options.prefix+'_read_delet.link')
			config.set('Read_link','chr_rr', options.prefix+'_read_chr_rr.link')
			config.set('Read_link','chr_rf', options.prefix+'_read_chr_rf.link')
			config.set('Read_link','chr_fr', options.prefix+'_read_chr_fr.link')
			config.set('Read_link','chr_ff', options.prefix+'_read_chr_ff.link')
			config.set('General','read_rf', 'yes')
			config.set('General','read_fr', 'yes')
			config.set('General','read_ff', 'yes')
			config.set('General','read_rr', 'yes')
			config.set('General','read_ins', 'yes')
			config.set('General','read_del', 'yes')
			config.set('General','read_chr_rr', 'yes')
			config.set('General','read_chr_rf', 'yes')
			config.set('General','read_chr_fr', 'yes')
			config.set('General','read_chr_ff', 'yes')
		else:
			config.set('General','read_rf', 'no')
			config.set('General','read_fr', 'no')
			config.set('General','read_ff', 'no')
			config.set('General','read_rr', 'no')
			config.set('General','read_ins', 'no')
			config.set('General','read_del', 'no')
			config.set('General','read_chr_rr', 'no')
			config.set('General','read_chr_rf', 'no')
			config.set('General','read_chr_fr', 'no')
			config.set('General','read_chr_ff', 'no')
			
		if options.dis_prop != 'not_filled':
			config.add_section('Proportion')
			config.set('Proportion','prop', options.dis_prop)
			config.set('General','prop', 'yes')
		else:
			config.set('General','prop', 'no')
			
		if options.agp != 'not_filled':
			config.add_section('Scaffold')
			config.set('Scaffold','scaff_tile', options.prefix+'_scaffold.tile')
			config.set('General','scaff_tile', 'yes')
		else:
			config.set('General','scaff_tile', 'no')
		# writting configuration file
		with open(options.prefix+'.conf', 'wb') as configfile:
			config.write(configfile)
	else:
		config.add_section('General')
		config.set('General','chr', options.chr)
		config.set('General','out_kar', options.out_kar)
		config.set('General','out_N', options.out_N)
		config.set('General','orient', options.orient)
		if options.cov != 'not_filled':
			config.add_section('Coverage')
			config.set('Coverage','cov', options.out_cov)
			config.set('General','cov', 'yes')
			config.set('Coverage','median_cov', str(Value[0]))
			config.set('Coverage','mean_cov', str(Value[1]))
		else:
			config.set('General','cov', 'no')
		config.add_section('Discord_link')
		config.add_section('Discord_zone')
		if options.frf != 'not_filled':
			config.set('Discord_link','frf', options.out_frf)
			config.set('General','frf', 'yes')
		else:
			config.set('General','frf', 'no')
			
		if options.ff != 'not_filled':
			config.set('Discord_link','ff', options.out_ff)
			config.set('General','ff', 'yes')
		else:
			config.set('General','ff', 'no')
			
		if options.rr != 'not_filled':
			config.set('Discord_link','rr', options.out_rr)
			config.set('General','rr', 'yes')
		else:
			config.set('General','rr', 'no')
			
		if options.ins != 'not_filled':
			config.set('Discord_link','ins', options.out_ins)
			config.set('General','ins', 'yes')
		else:
			config.set('General','ins', 'no')
			
		if options.delet != 'not_filled':
			config.set('Discord_link','delet', options.out_delet)
			config.set('General','delet', 'yes')
		else:
			config.set('General','delet', 'no')
			
		if options.chr_rr != 'not_filled':
			config.set('Discord_link','chr_rr', options.out_chr_rr)
			config.set('Discord_zone','chr_rr', options.chr_rr)
			config.set('General','chr_rr', 'yes')
		else:
			config.set('General','chr_rr', 'no')
			
		if options.chr_rf != 'not_filled':
			config.set('Discord_link','chr_rf', options.out_chr_rf)
			config.set('Discord_zone','chr_rf', options.chr_rf)
			config.set('General','chr_rf', 'yes')
		else:
			config.set('General','chr_rf', 'no')
			
		if options.chr_fr != 'not_filled':
			config.set('Discord_link','chr_fr', options.out_chr_fr)
			config.set('Discord_zone','chr_fr', options.chr_fr)
			config.set('General','chr_fr', 'yes')
		else:
			config.set('General','chr_fr', 'no')
			
		if options.chr_ff != 'not_filled':
			config.set('Discord_link','chr_ff', options.out_chr_ff)
			config.set('Discord_zone','chr_ff', options.chr_ff)
			config.set('General','chr_ff', 'yes')
		else:
			config.set('General','chr_ff', 'no')
			
		if options.liste_read != 'not_filled':
			config.add_section('Read_link')
			config.set('Read_link','rf', options.Rout_rf)
			config.set('Read_link','fr', options.Rout_fr)
			config.set('Read_link','ff', options.Rout_ff)
			config.set('Read_link','rr', options.Rout_rr)
			config.set('Read_link','ins', options.Rout_ins)
			config.set('Read_link','del', options.Rout_delet)
			config.set('Read_link','chr_rr', options.Rout_chr_rr)
			config.set('Read_link','chr_rf', options.Rout_chr_rf)
			config.set('Read_link','chr_fr', options.Rout_chr_fr)
			config.set('Read_link','chr_ff', options.Rout_chr_ff)
			config.set('General','read_rf', 'yes')
			config.set('General','read_fr', 'yes')
			config.set('General','read_ff', 'yes')
			config.set('General','read_rr', 'yes')
			config.set('General','read_ins', 'yes')
			config.set('General','read_del', 'yes')
			config.set('General','read_chr_rr', 'yes')
			config.set('General','read_chr_rf', 'yes')
			config.set('General','read_chr_fr', 'yes')
			config.set('General','read_chr_ff', 'yes')
		else:
			config.set('General','read_rf', 'no')
			config.set('General','read_fr', 'no')
			config.set('General','read_ff', 'no')
			config.set('General','read_rr', 'no')
			config.set('General','read_ins', 'no')
			config.set('General','read_del', 'no')
			config.set('General','read_chr_rr', 'no')
			config.set('General','read_chr_rf', 'no')
			config.set('General','read_chr_fr', 'no')
			config.set('General','read_chr_ff', 'no')
			
		if options.dis_prop != 'not_filled':
			config.add_section('Proportion')
			config.set('Proportion','prop', options.dis_prop)
			config.set('General','prop', 'yes')
		else:
			config.set('General','prop', 'no')
			
		if options.agp != 'not_filled':
			config.add_section('Scaffold')
			config.set('Scaffold','scaff_tile', options.out_scaff)
			config.set('General','scaff_tile', 'yes')
		else:
			config.set('General','scaff_tile', 'no')
		# writting configuration file
		with open(options.output, 'wb') as configfile:
			config.write(configfile)
		
if __name__ == "__main__": __main__()
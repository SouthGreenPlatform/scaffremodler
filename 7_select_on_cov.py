import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math, glob, datetime

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def run_job (cmd_line, ERROR):
	print cmd_line
	try:
		tmp = tempfile.NamedTemporaryFile().name
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

def run_job_silent (cmd_line, ERROR):
	# print cmd_line
	try:
		tmp = tempfile.NamedTemporaryFile().name
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


def mediane(L):
	L.sort()
	N = len(L)
	n = N/2.0
	p = int(n)
	if n == 0:
		return 0
	elif n == 1:
		return (L[0])
	elif n == p:
		return (L[p-1]+L[p])/2.0
	else:
		return float(L[p])
		
#########################################################################################################################################################
#				For calculating discordant coverage
#########################################################################################################################################################
def calcul_cov(LOCA_PROGRAMS, REF, SAM, TYPE, OUT):
	if TYPE == 'bam':
		cal_cov = '%s mpileup -A -f %s %s | cut -f 1,2,4 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), REF, SAM, OUT)
	elif TYPE == 'sam':
		cal_cov = '%s view -bS %s | %s mpileup -A -f %s - | cut -f 1,2,4 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, LOCA_PROGRAMS.get('Programs','samtools'), REF, OUT)
	else:
		mot = TYPE+' argument passed in --type is not recognized'
		sys.exit(mot)
	run_job(cal_cov, 'Error in calculating coverage:\n')

#########################################################################################################################################################
#				For zone identification
#########################################################################################################################################################


def select_sur_couv(COV, ZCOV, MAXCOV, MINCOV, NCOV, OUT):
	debut, av = 0, 0
	chr = ""
	liste = []
	outfile = open(OUT,'w')
	file = open(COV)
	i = 0
	for line in file:
		data = line.split()
		if data:
			POS = int(data[1])
			COV = int(data[2])
			if chr == '':#on initialise sur le premier chromosome
				chr = data[0]
				debut = POS
				av = POS
				liste.append(COV)
				DIC = set()
				DIC.add(POS)
			elif chr != data[0]: #on change de chromosome
				val_med = mediane(liste)
				if len(liste) >= int(ZCOV) and val_med <= float(MAXCOV) and val_med >= float(MINCOV):
					outfile.write('\t'.join([chr,str(debut),str(av),str(len(liste)),str(val_med)])+'\n')
					i += 1
				chr = data[0]
				debut = POS
				av = POS
				liste = []
				liste.append(COV)
				DIC = set()
				DIC.add(POS)
			elif (POS - av) > int(NCOV):#on change de zone, distance maximale pour une zone ou la couverture chute a 0
				val_med = mediane(liste)
				if len(liste) >= int(ZCOV) and val_med <= float(MAXCOV) and val_med >= float(MINCOV):
					outfile.write('\t'.join([chr,str(debut),str(av),str(len(liste)),str(val_med)])+'\n')
					i += 1
				debut = POS
				av = POS
				liste = []
				liste.append(COV)
				DIC = set()
				DIC.add(POS)
			else:
				liste.append(COV)
				av = POS
				DIC.add(POS)
	outfile.close()
	return i

#########################################################################################################################################################
#				For mate zone identification
#########################################################################################################################################################
def decoup_chr_bam(LOCA_PROGRAMS, BAM, CHR, OUT, TYPE):
	file = open(CHR)
	liste_to_remove = []
	for line in file:
		data = line.split()
		if data:
			debut = 1
			while debut <= int(data[1]):
				tempo1 = tempfile.NamedTemporaryFile()
				tempo2 = tempfile.NamedTemporaryFile()
				if (debut + 999999) < int(data[1]):
					tempo1.write(data[0]+'\t'+str(debut)+'\t'+str(debut + 999999)+'\n')
					fin = debut + 999999
				else:
					tempo1.write(data[0]+'\t'+str(debut)+'\t'+data[1]+'\n')
					fin = int(data[1])
				tempo1.flush()
				if TYPE == 'bam':
					bam2subbam = '%s view -L %s %s | cut -f 1 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), tempo1.name, BAM, tempo2.name)
				elif TYPE == 'sam':
					bam2subbam = '%s view -S -L %s %s | cut -f 1 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), tempo1.name, BAM, tempo2.name)
				else:
					mot = TYPE+' argument passed in --type is not recognized'
					sys.exit(mot)
				run_job_silent(bam2subbam, 'Error in decoup_chr_bam:\n')
				tempo2.flush()
				#Here we have read mapping in the zone but we need the second mate... that is why we use picar tools
				if os.path.getsize(tempo2.name) != 0:
					bam2subbam2 = '%s -jar %s INPUT=%s OUTPUT=%s QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool_FilterSamReads'), BAM, OUT+'_'+data[0]+'_'+str(debut)+'-'+str(fin)+'.bam', tempo2.name)
					run_job_silent(bam2subbam2, 'Error in decoup_chr_bam:\n')
				# print tempo1.name
				# print tempo2.name
				liste_to_remove.append(BAM+'_'+data[0]+'.bam')
				tempo1.close()
				tempo2.close()
				debut = debut + 500000
	return liste_to_remove

def create_sub_sam(LOCA_PROGRAMS, BAM, CHR, START, END, OUT):
	tempo1 = tempfile.NamedTemporaryFile()
	tempo1.write(CHR+'\t'+START+'\t'+END+'\n')
	tempo1.flush()
	tempo2 = tempfile.NamedTemporaryFile()
	bam2subbam = '%s view -bu -L %s %s | %s view -h -o %s -' % (LOCA_PROGRAMS.get('Programs','samtools'), tempo1.name, BAM, LOCA_PROGRAMS.get('Programs','samtools'), tempo2.name)
	run_job_silent(bam2subbam, 'Error in create_sub_sam:\n')
	tempo2.flush()
	if os.path.getsize(tempo2.name) != 0:
		bam2subbam2 = '%s -jar %s INPUT=%s OUTPUT=%s QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool_FilterSamReads'), BAM, OUT, tempo2.name)
		run_job_silent(bam2subbam2, 'Error in create_sub_sam:\n')
		# print tempo1.name
		# print tempo2.name
	else:
		print "warning there is a bug in create_sub_sam"
	tempo1.close()
	tempo2.close()

def trie_sam(SAM, OUT, CHR):
	outfile1 = open(OUT+'_tempo_amont','w')
	outfile2 = open(OUT+'_tempo_aval','w')
	outfile3 = open(OUT+'_tempo_chr','w')
	file = open(SAM)
	total = 0
	# amont = 0
	# aval = 0
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '@':
				outfile1.write(line)
				outfile2.write(line)
				outfile3.write(line)
			else:
				if data[2] != CHR:
					total += 1
					outfile3.write(line)
				elif int(data[3]) < int(data[7]):
					total += 1
					outfile1.write(line)
				elif int(data[3]) > int(data[7]):
					total += 1
					outfile2.write(line)
	outfile1.close()
	outfile2.close()
	outfile3.close()
	return total/2

def search_dest(LOCA_PROGRAMS, TARGET, DEST, DEST_CHR, ECART, START, END, REF, ZCOV, MAXCOV, MINCOV, OUT, TOTAL):
	outfile = open(OUT, 'a')
	#loading TARGET information
	file = open(TARGET)
	dico = set()
	for line in file:
		data = line.split()
		if data:
			if data[0][0] != '@':
				dico.add(line)
	file.close()
	while dico:
		size_liste = 0
		liste = []
		while size_liste == 0 and dico:
			line1 = dico.pop()
			decoupe = line1.split()
			if int(decoupe[3]) <= END and int(decoupe[3]) >= START:
				liste.append(line1)
				center = int(decoupe[7])
				L_dest = [center]
				chr_dest = decoupe[6]
				chr_target = decoupe[2]
				value = (((END - START) + 1) + ECART)
				deb_zone = (center - value)
				fin_zone = (center + value)
			size_liste = len(liste)
		size_liste = 0 #To initialise the loop
		while size_liste != len(liste):
			size_liste = len(liste)
			for n in dico:
				decoupe = n.split()
				if decoupe[6] == chr_dest and int(decoupe[3]) <= END and int(decoupe[3]) >= START:#The read is in the target zone and the destination is on the same chromosome
					if int(decoupe[7]) >= deb_zone and int(decoupe[7]) <= fin_zone:#The read has a similar destination
						L_dest.append(int(decoupe[7]))
						med_val = mediane(L_dest)
						value = (((END - START) + 1) + ECART)
						deb_zone = (med_val - value)
						fin_zone = (med_val + value)
						liste.append(n)
			for n in liste:
				if n in dico:
					dico.remove(n)
		if liste:
			tempo = tempfile.NamedTemporaryFile()
			for n in liste:
				data = n.split()
				tempo.write(data[0]+'\n')
			tempo.flush()
			info_target = recalc_border(LOCA_PROGRAMS, tempo.name, TARGET, REF, ZCOV, MAXCOV, MINCOV)
			if chr_dest == '=':
				info_dest = recalc_border(LOCA_PROGRAMS, tempo.name, DEST, REF, ZCOV, MAXCOV, MINCOV)
			else:
				info_dest = recalc_border(LOCA_PROGRAMS, tempo.name, DEST_CHR, REF, ZCOV, MAXCOV, MINCOV)
			tempo.close()
			if info_target[4] == 'PASS' and info_dest[4] == 'PASS':
				outfile.write('\t'.join([info_target[0], info_target[1], info_target[2], str(int(info_target[2])-int(info_target[1])+1), str(info_target[3]), info_dest[0], info_dest[1], info_dest[2], str(int(info_dest[2])-int(info_dest[1])+1), str(info_dest[3]), str(len(liste)/float(TOTAL)), str(len(liste))])+'\n')
	outfile.close()

def recalc_border(LOCA_PROGRAMS, LISTE, SAM, REF, ZCOV, MAXCOV, MINCOV):
	sub_sorting = '%s -jar %s INPUT=%s OUTPUT=%s QUIET=true FILTER=includeReadList READ_LIST_FILE=%s WRITE_READS_FILES=false VERBOSITY=WARNING VALIDATION_STRINGENCY=SILENT' % (LOCA_PROGRAMS.get('Programs','java'), LOCA_PROGRAMS.get('Programs','picard-tool_FilterSamReads'), SAM, SAM+'_sub.bam', LISTE)
	run_job_silent(sub_sorting, 'Error in recalc_border:\n')
	coverage = '%s mpileup -A -f %s %s | cut -f 1,2,4 > %s' % (LOCA_PROGRAMS.get('Programs','samtools'), REF, SAM+'_sub.bam', SAM+'_sub_cov')
	run_job_silent(coverage, 'Error in recalc_border:\n')
	chr = ""
	file = open(SAM+'_sub_cov')
	liste_cov = []
	for line in file:
		data = line.split()
		if data != []:
			liste_cov.append(int(data[2]))
			if chr == "":
				chr = data[0]
				min = data[1]
				max = data[1]
			elif chr == data[0]:
				max = data[1]
			else:
				sys.exit('There is a bug in recalc_border : several chromosomes are found in the cov file')
	MEDIAN = mediane(liste_cov)
	if len(liste_cov) >= ZCOV and MEDIAN >= MINCOV and MEDIAN <= MAXCOV:
		return [chr, min, max, MEDIAN, 'PASS']
	else:
		return [chr, min, max, MEDIAN, 'NO_PASS']


def look_4_mate(LOCA_PROGRAMS, TYPE, BAM, CHR, OUT, ZONE, ECART, REF, ZCOV, MAXCOV, MINCOV, NB_ZONE):
	#loading chromosome informations in a dictionnary
	dico_chr_info = {}
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			dico_chr_info[data[0]] = int(data[1])
	file.close()
	#for speed increase creation of a sub_bam file for each chromosome
	tmp = tempfile.NamedTemporaryFile().name
	current_dir_path = os.getcwd()
	out_tmp = current_dir_path+'/'+tmp.split('/')[-1]
	decoup_chr_bam(LOCA_PROGRAMS, BAM, CHR, out_tmp, TYPE)
	#Now it's time to work on each zone
	outfile = open(OUT,'w')
	outfile.close()
	file = open(ZONE)
	i = 0
	t0 = datetime.datetime.now()
	####This part to manage the sub_bam files
	chromosome = ""
	for line in file:
		data = line.split()
		if data != []:
			pos_zone_debut = int(data[1])
			pos_zone_fin = int(data[2])
			if chromosome != data[0]:
				chromosome = data[0]
				pos_debut = 1
				pos_fin = pos_debut + 999999
				if pos_fin > dico_chr_info[data[0]]:
					pos_fin = dico_chr_info[data[0]]
				while not((pos_debut <= pos_zone_debut) and (pos_zone_fin <= pos_fin)):
					pos_debut = pos_debut + 500000
					pos_fin = pos_debut + 999999
					if pos_fin > dico_chr_info[data[0]]:
						pos_fin = dico_chr_info[data[0]]
			else:
				while not((pos_debut <= pos_zone_debut) and (pos_zone_fin <= pos_fin)):
					pos_debut = pos_debut + 500000
					pos_fin = pos_debut + 999999
					if pos_fin > dico_chr_info[data[0]]:
						pos_fin = dico_chr_info[data[0]]
			#creation of a sub bam containing read mapping in each putative discordant zone and their mate
			create_sub_sam(LOCA_PROGRAMS, out_tmp+'_'+data[0]+'_'+str(pos_debut)+'-'+str(pos_fin)+'.bam', data[0], data[1], data[2], out_tmp+'_tempo_liste_'+data[0]+'.sam')
			#identifies mate mapping first and mate mapping second
			total = trie_sam(out_tmp+'_tempo_liste_'+data[0]+'.sam', out_tmp+'_.sam_'+data[0], data[0])
			#search for destination of second mate and perform selection
			search_dest(LOCA_PROGRAMS, out_tmp+'_.sam_'+data[0]+'_tempo_amont', out_tmp+'_.sam_'+data[0]+'_tempo_aval', out_tmp+'_.sam_'+data[0]+'_tempo_chr', float(ECART), int(data[1]), int(data[2]), REF, ZCOV, MAXCOV, MINCOV, OUT, total)
			search_dest(LOCA_PROGRAMS, out_tmp+'_.sam_'+data[0]+'_tempo_aval', out_tmp+'_.sam_'+data[0]+'_tempo_amont', out_tmp+'_.sam_'+data[0]+'_tempo_chr', float(ECART), int(data[1]), int(data[2]), REF, ZCOV, MAXCOV, MINCOV, OUT, total)
			i += 1
			if i % 100 == 0:
				print "Estimation for 7_step_remaining time for", BAM, ":", i ,"done in", datetime.datetime.now() - t0, ",",  NB_ZONE-i, "remaining"
	for filename in glob.glob(out_tmp+'_*'):
		# print filename
		os.remove(filename)

#########################################################################################################################################################
#				Merging identical and similar zones
#########################################################################################################################################################

def regroupe_filtre(DATA, DATA_prec, MAX):
	debut1 = str(min(int(DATA_prec[1]),int(DATA[1])))
	fin1 = str(max(int(DATA_prec[2]),int(DATA[2])))
	debut2 = str(min(int(DATA_prec[6]),int(DATA[6])))
	fin2 = str(max(int(DATA_prec[7]),int(DATA[7])))
	size1 = str((int(fin1)-int(debut1)) + 1)
	size2 = str((int(fin2)-int(debut2)) + 1)
	mean_cov1 = ((int(DATA_prec[12])*float(DATA_prec[4]))+(int(DATA[12])*float(DATA[4])))/(int(DATA_prec[12])+int(DATA[12]))
	mean_cov2 = ((int(DATA_prec[13])*float(DATA_prec[9]))+(int(DATA[13])*float(DATA[9])))/(int(DATA_prec[13])+int(DATA[13]))
	if (DATA[0] == DATA_prec[0]) and (DATA[5] == DATA_prec[5]):
		if int(DATA_prec[1]) <= int(DATA[1]):
			if ((int(DATA[1]) - int(DATA_prec[2])) < MAX) and ((int(DATA[1]) < int(DATA[6]) and int(DATA_prec[1]) < int(DATA_prec[6])) or (int(DATA[1]) > int(DATA[6]) and int(DATA_prec[1]) > int(DATA_prec[6]))):
				if int(DATA[7]) <= int(DATA_prec[6]):#la zone arrive avant la zone precedente : ssssss pppppp
					if ((int(DATA_prec[6])-int(DATA[7])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA_prec[7]) <= int(DATA[6]):#la zone arrive apres la zone precedente : pppppp ssssss
					if ((int(DATA[6])-int(DATA_prec[7])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA[7]) >= int(DATA_prec[6]) and int(DATA[7]) <= int(DATA_prec[7]) and int(DATA[6]) <= int(DATA_prec[6]): #zone precedente chevauche sur la suivante en mappant apres : sssspspspppppp
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA[7]) >= int(DATA_prec[7]) and int(DATA[6]) <= int(DATA_prec[6]): #la zone precedente est inclue dans la suivante  sssspspspspssss
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA_prec[7]) >= int(DATA[6]) and int(DATA_prec[7]) <= int(DATA[7]) and int(DATA[6]) >= int(DATA_prec[6]): #zone precedente chevauche sur la suivante en mappant avant : pppppspspspssssss
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA_prec[7]) >= int(DATA[7]) and int(DATA_prec[6]) <= int(DATA[6]): #la zone suivante est inclue dans la suivante  ppppppspspspspppppp
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				else:
					return [DATA_prec, 'not_found']
			else:
				return [DATA_prec, 'not_found']
		else: #DATA mapp before DATA_prec
			if ((int(DATA_prec[1]) - int(DATA[2])) < MAX) and ((int(DATA[1]) < int(DATA[6]) and int(DATA_prec[1]) < int(DATA_prec[6])) or (int(DATA[1]) > int(DATA[6]) and int(DATA_prec[1]) > int(DATA_prec[6]))):
				if int(DATA[7]) <= int(DATA_prec[6]):#la zone arrive avant la zone precedente : ssssss pppppp
					if ((int(DATA_prec[6])-int(DATA[7])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA_prec[7]) <= int(DATA[6]):#la zone arrive apres la zone precedente : pppppp ssssss
					if ((int(DATA[6])-int(DATA_prec[7])) < MAX):#on est sur deux zones contigues non chevauchantes
						return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
					else:#les zones ne sont pas contigues
						return [DATA_prec, 'not_found']
				elif int(DATA[7]) >= int(DATA_prec[6]) and int(DATA[7]) <= int(DATA_prec[7]) and int(DATA[6]) <= int(DATA_prec[6]): #zone precedente chevauche sur la suivante en mappant apres : sssspspspppppp
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA[7]) >= int(DATA_prec[7]) and int(DATA[6]) <= int(DATA_prec[6]): #la zone precedente est inclue dans la suivante  sssspspspspssss
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA_prec[7]) >= int(DATA[6]) and int(DATA_prec[7]) <= int(DATA[7]) and int(DATA[6]) >= int(DATA_prec[6]): #zone precedente chevauche sur la suivante en mappant avant : pppppspspspssssss
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				elif int(DATA_prec[7]) >= int(DATA[7]) and int(DATA_prec[6]) <= int(DATA[6]): #la zone suivante est inclue dans la suivante  ppppppspspspspppppp
					return [[DATA[0], debut1, fin1, size1, mean_cov1, DATA_prec[5], debut2, fin2, size2, mean_cov2, ((float(DATA_prec[10])*int(DATA_prec[11]))+(float(DATA[10])*int(DATA[11])))/(int(DATA_prec[11])+int(DATA[11])), int(DATA_prec[11])+int(DATA[11]), int(DATA_prec[12])+int(DATA_prec[3]), int(DATA_prec[13])+int(DATA_prec[8]), DATA[14]], 'found']
				else:
					return [DATA_prec, 'not_found']
			else:
				return [DATA_prec, 'not_found']
	else:
		return [DATA_prec, 'not_found']

def merge_zone(FILE, MAX, OUT):
	#recording each zone in a dictionary
	data_prec = ''
	file = open(FILE)
	outfile = open(OUT,'w')
	dic = set()
	i = 0
	len_dic = 0
	for line in file:
		data = line.split()
		if data:
			i = i + 1
			if not(i in dic):
				data_prec = line.split()
				data_prec.append(data_prec[3])
				data_prec.append(data_prec[8])
				data_prec.append(i)
				dic.add(i)
				while len(dic) != len_dic:
					len_dic = len(dic)
					file2 = open(FILE)
					j = 0
					for line2 in file2:
						data_check = line2.split()
						if data_check:
							j = j + 1
							if not(j in dic):
								data2 = line2.split()
								data2.append(data_prec[3])
								data2.append(data_prec[8])
								data2.append(j)
								found = regroupe_filtre(data2, data_prec, MAX)
								data_prec = found[0]
								#We try the other order of zone 'amont and aval exchange places'
								if found[1] == 'not_found':
									new_data2 = [data2[5],data2[6],data2[7],data2[8],data2[9],data2[0],data2[1],data2[2],data2[3],data2[4],data2[10],data2[11],data2[13],data2[12],data2[14]]
									found = regroupe_filtre(new_data2, data_prec, MAX)
									data_prec = found[0]
								elif found[1] != 'found':
									sys.exit('bug')
								dic.add(data_prec[14])
				outfile.write(str(data_prec[0])+' '+str(data_prec[1])+' '+str(data_prec[2])+' '+str(data_prec[3])+' '+str(data_prec[4])+' '+str(data_prec[5])+' '+str(data_prec[6])+' '+str(data_prec[7])+' '+str(data_prec[8])+' '+str(data_prec[9])+' '+str(data_prec[10])+' '+str(data_prec[11])+'\n')
	outfile.close()

#########################################################################################################################################################
#				Calculating score
#########################################################################################################################################################
def calculate_score(FILE, YIS, MIS, YIC, MIC, MIN, CHR, OUT):
	############################################
	#recording chromosome order
	############################################
	liste_chr = []
	file = open(CHR)
	for line in file:
		data = line.split()
		if data:
			liste_chr.append(data[0])
	file.close()
	
	#parametres affine sur multiplicateur de taille
	b = float(YIS)
	x = float(MIS)
	a = (1-b)/x

	#parametres affine sur multiplicateur de couverture
	d = float(YIC)
	y = float(MIC)
	c = (1-d)/y
	
	outfile = open(OUT, 'w')
	file = open(FILE)
	for line in file:
		data = line.split()
		if data:
			#score de depart
			score = 100
			#multiplicateur du score sur la taille
			if ((int(data[3])+int(data[8]))/2.0) >= x:
				score = score*1
			else:
				score = score*(a*((int(data[3])+int(data[8]))/2.0)+b)
			#multiplicateur du score sur la couverture
			if ((float(data[4])+float(data[9]))/2.0) >= y:
				score = score*1
			else:
				score = score*(c*((float(data[4])+float(data[9]))/2.0)+d)
			###############
			#Ordering data#
			###############
			if data[0] == data[5]:
				if int(data[1]) <= int(data[6]):
					data_out = list(data)
				else:
					data_out = [data[5], data[6], data[7], data[8], data[9], data[0], data[1], data[2], data[3], data[4], data[10], data[11]]
			else:
				if liste_chr.index(data[0]) < liste_chr.index(data[5]):
					data_out = list(data)
				else:
					data_out = [data[5], data[6], data[7], data[8], data[9], data[0], data[1], data[2], data[3], data[4], data[10], data[11]]
			data_out.append(str(score))
			if score >= MIN:
				data_out.append('PASSED')
			else:
				data_out.append('NOT_PASSED')
			outfile.write('\t'.join(data_out)+'\n')
	outfile.close()

	
def __main__():
	t_start = datetime.datetime.now()
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\nThis program take in input a sam/bam file,"
	" and calculate coverage for covered sites. The bam should be coordinate sorted")
	# Wrapper options. 
	parser.add_option( '', '--ref', dest='ref', default='not_filled', help='The multifasta reference file')
	parser.add_option( '', '--sam', dest='sam', default='not_filled', help='Paired sam/bam file')
	parser.add_option( '', '--type', dest='type', default='bam', help='Input type : sam or bam, [default: %default]')
	parser.add_option( '', '--min_zone', dest='min_zone', default=500, help='Minimal number of covered sites in the zone, [default: %default]')
	parser.add_option( '', '--maxcov', dest='maxcov', default=300, help='The maximal median coverage accepted (float), [default: %default]')
	parser.add_option( '', '--mincov', dest='mincov', default=0, help='The minimal median coverage accepted (float), [default: %default]')
	parser.add_option( '', '--min_gap', dest='min_gap', default=300, help='The maximal distance of contiguous uncovered sites, [default: %default]')
	parser.add_option( '', '--ecart', dest='ecart', default=2000, help='The value that will be added or substracted to keep a read as the same destination. Recommended 3*sd(insert)')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='File containing col1: chromosome name and col2: chromsome size')
	parser.add_option( '', '--max_dist_merge', dest='max_dist_merge', default=1000, help='Maximal distance between two discordant zone to merge, [default: %default]')
	parser.add_option( '', '--YiS', dest='YiS', default=0, help='The Y-intercept of the linear function for zone size that will give the first component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiS', dest='MiS', default=1000, help='The minimal zone size for which the first component of product giving the score will be maximal (integer), [default: %default]')
	parser.add_option( '', '--YiC', dest='YiC', default=0, help='The Y-intercept of the linear function for coverage that will give the second component of product giving the score (integer), [default: %default]')
	parser.add_option( '', '--MiC', dest='MiC', default=25, help='The minimal zone coverage for which the second component of product giving the score will be maximal (integer), [default: %default]')
	parser.add_option( '', '--min_score', dest='min_score', default=70, help='The minimal score for a discordant zone to be identified as passed, [default: %default]')
	parser.add_option( '', '--out', dest='out', default='not_filled', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	pathname = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	
	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		sys.exit(mot)
	if options.sam == 'not_filled':
		mot = 'Please provide an argument for --sam'
		sys.exit(mot)
	
	if os.path.getsize(options.sam) == 0:
		outfile = open(options.out, 'w')
		outfile.close()
	else:
		tmp_name = tempfile.NamedTemporaryFile().name
		tmp_cov = tmp_name+'.cov'
		tmp_zone = tmp_name+'.zone'
		tmp_mate_zone = tmp_name+'_mate.zone'
		tmp_merge = tmp_name+'.merge'
		if options.config:
			config = ConfigParser.RawConfigParser()
			config.read(options.config)
			calcul_cov(loca_programs, config.get('General','ref'), options.sam, config.get('Trie_discord','type'), tmp_cov)
			maxcov = config.getfloat('Calc_coverage','median_coverage')*config.getfloat('General','mult_max_cov')
			mincov = config.getfloat('Calc_coverage','median_coverage')*config.getfloat('General','mult_min_cov')
			print 'Minimal accepted coverage:', mincov
			print 'Maximal accepted coverage:', maxcov
			nb_zone = select_sur_couv(tmp_cov, config.getint('General','min_zone'), maxcov, mincov, config.getint('General','min_gap'), tmp_zone)
			os.remove(tmp_cov)
			ecart = config.getfloat('Calc_coverage','standard_deviation_insert')*3.0
			print 'Margin:', ecart
			print 'Number of zone to test:', nb_zone
			look_4_mate(loca_programs, config.get('Trie_discord','type'), options.sam, config.get('General','chr'), tmp_mate_zone, tmp_zone, ecart, config.get('General','ref'), config.getint('General','min_zone'), maxcov, mincov, nb_zone)
			os.remove(tmp_zone)
			merge_zone(tmp_mate_zone, config.getfloat('General','max_dist_merge'), tmp_merge)
			os.remove(tmp_mate_zone)
			calculate_score(tmp_merge, config.getfloat('General','YiS'), config.getfloat('Score_discord','MiS'), config.getfloat('General','YiC'), config.getfloat('Score_discord','MiC'), config.getfloat('General','min_score'), config.get('General','chr'), options.out)
			os.remove(tmp_merge)
			os.system("sed -i '1i#CHR-zone1\tSTART\tEND\tSIZE\tCOV\tCHR-zone2\tSTART\tEND\SIZE\tCOV\tMISC\tREAD\tSCORE\tSTATUS' %s" % options.out)
		else:
			if options.ref == 'not_filled':
				mot = 'Please provide an argument for --ref'
				sys.exit(mot)
			calcul_cov(loca_programs, options.ref, options.sam, options.type, tmp_cov)
			nb_zone = select_sur_couv(tmp_cov, int(options.min_zone), float(options.maxcov), float(options.mincov), int(options.min_gap), tmp_zone)
			print 'Number of zone to test:', nb_zone
			os.remove(tmp_cov)
			look_4_mate(loca_programs, options.type, options.sam, options.chr, tmp_mate_zone, tmp_zone, float(options.ecart), options.ref, int(options.min_zone), float(options.maxcov), float(options.mincov), nb_zone)
			os.remove(tmp_zone)
			merge_zone(tmp_mate_zone, int(options.max_dist_merge), tmp_merge)
			os.remove(tmp_mate_zone)
			calculate_score(tmp_merge, float(options.YiS), float(options.MiS), float(options.YiC), float(options.MiC), float(options.min_score), options.chr, options.out)
			os.remove(tmp_merge)
			os.system("sed -i '1i#CHR-zone1\tSTART\tEND\tSIZE\tCOV\t#CHR-zone2\tSTART\tEND\SIZE\tCOV\tMISC\tREAD\tSCORE\tSTATUS' %s" % options.out)
	print "Temps total :", datetime.datetime.now() - t_start
if __name__ == "__main__": __main__()
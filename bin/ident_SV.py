import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, random

# def couv2couv(CHR, COV, PREFIX):
	# liste = []
	# file = open(CHR)
	# for line in file:
		# data = line.split()
		# if data:
			# grep(COV, PREFIX+'_'+data[0], data[0], 'keep', 'all')
			# liste.append(PREFIX+'_'+data[0])
	# return liste

# def grep(FILE, OUT, MOT, KEEP, COL):
	# outfile = open(OUT, 'w')
	# i = 0
	# if COL == 'all':
		# if KEEP == 'keep':
			# for line in open(FILE):
				# if MOT in line.split():
					# outfile.write(line)
					# i += 1
		# elif KEEP == 'rm':
			# for line in open(FILE):
				# if not(MOT in line.split()):
					# outfile.write(line)
					# i += 1
		# else:
			# sys.exit('bug in grep')
	# else:
		# if KEEP == 'keep':
			# for line in open(FILE):
				# if MOT in line.split():
					# data = line.split()
					# liste = []
					# for n in COL:
						# liste.append(data[n])
					# outfile.write('\t'.join(liste)+'\n')
					# i += 1
		# elif KEEP == 'rm':
			# for line in open(FILE):
				# if not(MOT in line.split()):
					# data = line.split()
					# liste = []
					# for n in COL:
						# liste.append(data[n])
					# outfile.write('\t'.join(liste)+'\n')
					# i += 1
		# else:
			# sys.exit('bug in grep')
	# outfile.close()
	# return i
	
def couv2couv(chrFile, covFile, covDic):
	"""
		Stock in memory the coverture of each chromosomes
		
		:param chrFile: the path to the tabular chromosome file : chr_name chr_size
		:type chrFile: str:
		:param chrFile: the path to the tabular coverture file : chr_name position coverture
		:type chrFile: str:
		:param covDic: The dictionnary to stock the coverture
		:type covDic: dict
	"""
	chr = ""
	chrSize = {}
	tempListeCov = []
	i = 0
	chrInput = open(chrFile, 'r')
	for line in chrInput:  # Stock the size of each chromosome in a dictionnary
		if line.strip():
			cols = line.split()
			chrSize[cols[0]] = int(cols[1])
	chrInput.close()
	
	covInput = open(covFile, 'r')
	for line in covInput:
		if line.strip():
			cols = line.split()
			if not cols[0] in covDic:  # new chromosome
				if not chr:  # First line
					while i < int(cols[1]):
						tempListeCov.append(0)
						i += 1
					tempListeCov.append(int(cols[2]))
					i += 1
					chr = cols[0]
				else:
					while i < chrSize[chr]:
						tempListeCov.append(0)
						i += 1
					covDic[chr] = tempListeCov
					tempListeCov = []
					chr = cols[0]
					i = 0
					while i < int(cols[1]):
						tempListeCov.append(0)
						i += 1
					tempListeCov.append(int(cols[2]))
					i += 1
			else:  # same chromosome
				while i < int(cols[1]):
					tempListeCov.append(0)
					i += 1
				tempListeCov.append(int(cols[2]))
				i += 1
	while i < chrSize[chr]:
		tempListeCov.append(0)
		i += 1
	covDic[chr] = tempListeCov
	tempListeCov = []

def indent_discord(FF, FR, RR, INS, DEL, CHR_rr, CHR_fr, CHR_rf, CHR_ff, INSERT, OUT, EXP_COV, PLOID, TYPE):
	outfile = open(OUT,'w')
	file = open(FF)
	dic_FF = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_FF:
					dic_FF[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_FF[data[0]] = ([])
					dic_FF[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(FR)
	dic_FR = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_FR:
					dic_FR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_FR[data[0]] = ([])
					dic_FR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(RR)
	dic_RR = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_RR:
					dic_RR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_RR[data[0]] = ([])
					dic_RR[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(INS)
	dic_INS = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_INS:
					dic_INS[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_INS[data[0]] = ([])
					dic_INS[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(DEL)
	dic_DEL = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_DEL:
					dic_DEL[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_DEL[data[0]] = ([])
					dic_DEL[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
	file.close()
	
	file = open(CHR_rr)
	dic_CHR_rr = {}
	dic_CHR_rr_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_rr:
					dic_CHR_rr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_rr[data[0]] = ([])
					dic_CHR_rr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_rr_rev:
					dic_CHR_rr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_rr_rev[data[5]] = ([])
					dic_CHR_rr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	
	file = open(CHR_rf)
	dic_CHR_rf = {}
	dic_CHR_rf_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_rf:
					dic_CHR_rf[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_rf[data[0]] = ([])
					dic_CHR_rf[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_rf_rev:
					dic_CHR_rf_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_rf_rev[data[5]] = ([])
					dic_CHR_rf_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	file = open(CHR_fr)
	dic_CHR_fr = {}
	dic_CHR_fr_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_fr:
					dic_CHR_fr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_fr[data[0]] = ([])
					dic_CHR_fr[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_fr_rev:
					dic_CHR_fr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_fr_rev[data[5]] = ([])
					dic_CHR_fr_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	
	file = open(CHR_ff)
	dic_CHR_ff = {}
	dic_CHR_ff_rev = {}
	for line in file:
		data = line.split()
		if data != []:
			if data[0][0] != "#" and data[13] == "PASSED":
				if data[0] in dic_CHR_ff:
					dic_CHR_ff[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				else :
					dic_CHR_ff[data[0]] = ([])
					dic_CHR_ff[data[0]].append([int(data[1]), int(data[2]), data[5], int(data[6]), int(data[7])])
				if data[5] in dic_CHR_ff_rev:
					dic_CHR_ff_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
				else :
					dic_CHR_ff_rev[data[5]] = ([])
					dic_CHR_ff_rev[data[5]].append([int(data[6]), int(data[7]), data[0], int(data[1]), int(data[2])])
	file.close()
	vasouille = 100
	#All discordant position has been loaded
	##################################################################
	#####This part only works for intra chromosomal rearrangments#####
	##################################################################
	#Now it's time to detect simple translocations
	if TYPE == '0':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if k[3] <= pos3[1] + 2*INSERT and  k[3] + vasouille >= pos3[1] and k[0] + vasouille >= pos1[1] and k[1] <= pos3[0] + vasouille:#potential translocation
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for m in dic_DEL[n]:
								if m[1] <= pos1[0] + vasouille and m[1] + 2*INSERT >= pos1[0] and m[3] + vasouille >= pos2[1] and m[3] <= pos2[1] + 2*INSERT:
									pos0 = [m[0],m[1]]
									pos5 = [m[3],m[4]]
									outfile.write('translocation\tregion:\t'+n+'\t'+str(pos5[0])+'\t'+str(pos3[1])+'\ttarget:\t'+n+'\t'+str(pos0[1])+'\t'+str(pos1[0])+'\tALTERNATIVE\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos3[1])+'\t'+str(pos4[0])+'\n')
	#Now it's time to detect simple duplications
	if TYPE == '1':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if k[3] <= pos3[1] + 2*INSERT and  k[3] + vasouille >= pos3[1] and k[0] + vasouille >= pos1[1] and k[1] <= pos3[0] + vasouille:#potential duplication
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							#Maintenant s'occuper des couvertures
							if couv_med_reg(pos1[0], pos2[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos3[1])+'\t'+str(pos4[0])+'\n')
						elif pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] + vasouille >= k[4] and pos1[1] <= k[3] + vasouille:#potential duplication
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							#Maintenant s'occuper des couvertures
							if couv_med_reg(pos4[0], pos3[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+n+'\t'+str(pos4[0])+'\t'+str(pos3[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect simple inversions
	if TYPE == '0':
		for n in dic_RR:
			for j in dic_RR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_FF:
					for k in dic_FF[n]:
						if k[0] <= pos1[1] + 2*INSERT and  k[0] + vasouille >= pos1[1] and k[3] <= pos3[1] + 2*INSERT and k[3] + vasouille >= pos3[1]:#simple inversion
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							outfile.write('inversion\tregion:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos3[1])+'\n')
	#Now it's time to detect simple deletions
	if TYPE == '2':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if couv_med_reg_del(pos1[1], pos3[0], n) <= (EXP_COV - float(EXP_COV)*PLOID):
					outfile.write('deletion\tregion:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple tandem duplication
	if TYPE == '3':
		for n in dic_FR:
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if couv_med_reg(pos1[0], pos3[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
					outfile.write('tandem_duplication\tregion:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos3[1])+'\n')
	#Now it's time to detect simple reciprocal translocations
	if TYPE == '0':
		for n in dic_FR:
			dic_DEL_FR = []
			for j in dic_FR[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				if n in dic_DEL:
					for k in dic_DEL[n]:
						if pos1[0] <= k[1] + 2*INSERT and pos1[0] + vasouille >= k[1] and k[3] <= pos3[1] + 2*INSERT and k[3] + vasouille >= pos3[1]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							dic_DEL_FR.append([pos2[0], pos1[1], pos3[0], pos4[1], pos1[0], pos4[0], pos2[1], pos3[1]])
			L_done = []
			for x in dic_DEL_FR:
				L_done.append(x)
				for y in dic_DEL_FR:
					if not(y in L_done):
						if x[1] <= y[0] + vasouille and y[1] <= x[2] + vasouille and x[3] <= y[2] + vasouille:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(x[4])+'\t'+str(y[6])+'\tregion2:\t'+n+'\t'+str(x[5])+'\t'+str(y[7])+'\n')
						elif y[1] <= x[0] + vasouille and x[1] <= y[2] + vasouille and y[3] <= x[2] + vasouille:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(y[4])+'\t'+str(x[6])+'\tregion2:\t'+n+'\t'+str(y[5])+'\t'+str(x[7])+'\n')
	#Now it's time to detect inversion of a translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[0] + vasouille >= pos1[1]:# F RR F structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_DEL:
								for m in dic_DEL[n]:
									if m[1] <= pos1[0] + vasouille and m[1] + 2*INSERT >= pos1[0] and m[3] + vasouille >= pos2[1] and m[3] <= pos2[1] + 2*INSERT:# DEL F R DEL R F structure
										pos0 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos4[1])+'\t'+str(pos3[0])+'\n')
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and k[3] + vasouille >= pos3[1]:# R FF R structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_DEL:
								for m in dic_DEL[n]:
									if m[1] <= pos3[0] + vasouille and m[1] + 2*INSERT >= pos3[0] and m[3] + vasouille >= pos4[1] and m[3] <= pos4[1] + 2*INSERT:# R F DEL F R DEL structure
										pos0 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect inversion of a duplication
	if TYPE == '4':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[0] + vasouille >= pos1[1]:# F RR F structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos1[0], pos2[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos4[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '5':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and k[3] + vasouille >= pos3[1]:# R FF R structure
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos3[0], pos4[1], n) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos2[1])+'\t'+str(pos1[0])+'\n')
	#Now it's time to detect inversion of both fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_RR:
					for k in dic_RR[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for x in dic_FF[n]:
								if x != j:
									pos5 = [x[0], x[1]]
									pos6 = [x[3], x[4]]
									for y in dic_RR[n]:
										if y != k:
											if pos5[0] <= y[1] + 2*INSERT and  pos5[0] + vasouille >= y[1] and pos6[0] <= y[4] + 2*INSERT and  pos6[0] + vasouille >= y[4]:
												if pos1[1] <= y[0] + vasouille and pos6[1] <= pos4[0] + vasouille:
													outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(y[1])+'\tregion2_inv:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\n')
	#Now it's time to detect inversion of the first or the second fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_FF:
			for j in dic_FF[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_FR:
					for k in dic_FR[n]:
						if pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_RR:
								for x in dic_RR[n]:
									if pos2[0] + vasouille >= x[1] and pos2[0] <= x[1] + 2*INSERT:
										pos5 = [x[0], x[1]]
										pos6 = [x[3], x[4]]
										if n in dic_DEL:
											for y in dic_DEL[n]:
												if y[3] + vasouille >= pos6[1] and y[3] <= pos6[1] + 2*INSERT and y[1] <= pos1[0] + vasouille and y[1] + 2*INSERT >= pos1[0]:
													pos7 = [y[0], y[1]]
													pos8 = [y[3], y[4]]
													if pos2[0] + vasouille >= pos1[1]:#first fragment inversed
														outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos5[1])+'\tregion2:\t'+n+'\t'+str(pos8[0])+'\t'+str(pos4[1])+'\n')
													elif pos2[1] <= pos1[0] + vasouille:#second fragment inversed
														outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos7[1])+'\tregion2_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos6[1])+'\n')
	##################################################################
	#####This part only works for inter chromosomal rearrangments#####
	##################################################################
	#Now it's time to detect simple translocations
	if TYPE == '0':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				pos5 = []
				pos6 = []
				if n in dic_CHR_fr_rev:
					for k in dic_CHR_fr_rev[n]:
						if k[1] <= pos2[0] + vasouille and k[1] + 2*INSERT >= pos2[0]:
							pos4 = [k[0], k[1]]
							pos3 = [k[3], k[4]]
							if n in dic_CHR_rf_rev:
								for m in dic_CHR_rf_rev[n]:
									if m[0] + vasouille >= pos1[1] and m[0] <= pos1[1] + 2*INSERT and m[2] == k[2]:
										pos6 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										if pos5[1] <= pos3[0] + vasouille and pos5[1] + 2*INSERT >= pos3[0]:
											outfile.write('translocation\tregion:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\ttarget:\t'+k[2]+'\t'+str(pos5[1])+'\t'+str(pos3[0])+'\n')
				if n in dic_CHR_fr:
					for k in dic_CHR_fr[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT:
							pos4 = [k[0], k[1]]#inversion pos4 et pos3
							pos3 = [k[3], k[4]]
							if n in dic_CHR_rf:
								for m in dic_CHR_rf[n]:
									if m[1] <= pos2[0] + vasouille and m[1] + 2*INSERT >= pos2[0] and m[2] == k[2]:
										pos6 = [m[0],m[1]]#inversion pos6 et pos5
										pos5 = [m[3],m[4]]
										if pos3[1] <= pos5[0] + vasouille and pos3[1] + 2*INSERT >= pos5[0]:
											outfile.write('translocation\tregion:\t'+n+'\t'+str(pos4[0])+'\t'+str(pos6[1])+'\ttarget:\t'+k[2]+'\t'+str(pos3[1])+'\t'+str(pos5[0])+'\n')
	#Now it's time to detect simple duplications
	if TYPE == '6':
		for n in dic_CHR_rf:
			for j in dic_CHR_rf[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_fr:
					for k in dic_CHR_fr[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos2[0], pos4[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+k[2]+'\t'+str(pos2[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '7':
		for n in dic_CHR_fr_rev:
			for j in dic_CHR_fr_rev[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_rf_rev:
					for k in dic_CHR_rf_rev[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos2[0], pos4[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion:\t'+k[2]+'\t'+str(pos2[0])+'\t'+str(pos4[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple translocations with inversion
	if TYPE == '0':
		for n in dic_DEL:
			for j in dic_DEL[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				pos5 = []
				pos6 = []
				if n in dic_CHR_rr_rev:
					for k in dic_CHR_rr_rev[n]:
						if k[1] <= pos2[0] + vasouille and k[1] + 2*INSERT >= pos2[0]:
							pos4 = [k[0], k[1]]
							pos3 = [k[3], k[4]]
							if n in dic_CHR_ff_rev:
								for m in dic_CHR_ff_rev[n]:
									if m[0] + vasouille >= pos1[1] and m[0] <= pos1[1] + 2*INSERT and m[2] == k[2]:
										pos6 = [m[0],m[1]]
										pos5 = [m[3],m[4]]
										if pos3[1] <= pos5[0] + vasouille and pos3[1] + 2*INSERT >= pos5[0]:
											outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\ttarget:\t'+k[2]+'\t'+str(pos3[1])+'\t'+str(pos5[0])+'\n')
				if n in dic_CHR_ff:
					for k in dic_CHR_ff[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if n in dic_CHR_rr:
								for m in dic_CHR_rr[n]:
									if m[1] <= pos2[0] + vasouille and m[1] + 2*INSERT >= pos2[0] and m[2] == k[2]:
										pos5 = [m[0],m[1]]
										pos6 = [m[3],m[4]]
										if pos6[1] <= pos4[0] + vasouille and pos6[1] + 2*INSERT >= pos4[0]:
											outfile.write('translocation\tregion_inv:\t'+n+'\t'+str(pos3[0])+'\t'+str(pos5[1])+'\ttarget:\t'+k[2]+'\t'+str(pos6[1])+'\t'+str(pos4[0])+'\n')
	#Now it's time to detect simple duplications with inversion
	if TYPE == '8':
		for n in dic_CHR_rr:
			for j in dic_CHR_rr[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_ff:
					for k in dic_CHR_ff[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos4[0], pos2[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+k[2]+'\t'+str(pos4[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	if TYPE == '9':
		for n in dic_CHR_rr_rev:
			for j in dic_CHR_rr_rev[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				pos3 = []
				pos4 = []
				if n in dic_CHR_ff_rev:
					for k in dic_CHR_ff_rev[n]:
						if k[0] + vasouille >= pos1[1] and k[0] <= pos1[1] + 2*INSERT and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							if couv_med_reg(pos4[0], pos2[1], j[2]) >= (EXP_COV + float(EXP_COV)*PLOID):
								outfile.write('duplication\tregion_inv:\t'+k[2]+'\t'+str(pos4[0])+'\t'+str(pos2[1])+'\ttarget:\t'+n+'\t'+str(pos1[1])+'\t'+str(pos3[0])+'\n')
	#Now it's time to detect simple reciprocal translocations
	if TYPE == '0':
		for n in dic_CHR_fr:
			dic_CHR_fr_rf = []
			for j in dic_CHR_fr[n]:
				pos1 = [j[0], j[1]]
				pos2 = [j[3], j[4]]
				if n in dic_CHR_rf:
					for k in dic_CHR_rf[n]:
						if pos1[0] <= k[1] + 2*INSERT and pos1[0] + vasouille >= k[1] and k[3] <= pos2[1] + 2*INSERT and k[3] + vasouille >= pos2[1] and j[2] == k[2]:
							pos3 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							dic_CHR_fr_rf.append([pos3[0], pos1[1], pos2[0], pos4[1], pos1[0], pos4[0], pos3[1], pos2[1], j[2]])
			L_done = []
			for x in dic_CHR_fr_rf:
				L_done.append(x)
				for y in dic_CHR_fr_rf:
					if not(y in L_done):
						if x[1] <= y[0] + vasouille and  x[3] <= y[2] + vasouille and x[8] == y[8]:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(x[4])+'\t'+str(y[6])+'\tregion2:\t'+x[8]+'\t'+str(x[5])+'\t'+str(y[7])+'\n')
						elif y[1] <= x[0] + vasouille and y[3] <= x[2] + vasouille and x[8] == y[8]:
							outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(y[4])+'\t'+str(x[6])+'\tregion2:\t'+x[8]+'\t'+str(y[5])+'\t'+str(x[7])+'\n')
	#Now it's time to detect inversion of both fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_CHR_ff:
			for j in dic_CHR_ff[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if n in dic_CHR_rr:
					for k in dic_CHR_rr[n]:
						if pos1[0] <= k[1] + 2*INSERT and  pos1[0] + vasouille >= k[1] and pos3[0] <= k[4] + 2*INSERT and  pos3[0] + vasouille >= k[4] and k[2] == j[2]:
							pos2 = [k[0], k[1]]
							pos4 = [k[3], k[4]]
							for x in dic_CHR_ff[n]:
								if x != j:
									pos5 = [x[0], x[1]]
									pos6 = [x[3], x[4]]
									for y in dic_CHR_rr[n]:
										if y != k:
											if pos5[0] <= y[1] + 2*INSERT and  pos5[0] + vasouille >= y[1] and pos6[0] <= y[4] + 2*INSERT and  pos6[0] + vasouille >= y[4] and x[2] == y[2] and x[2] == k[2]:
												if pos1[1] <= y[0] + vasouille and pos6[1] <= pos4[0] + vasouille:
													outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(y[1])+'\tregion2_inv:\t'+k[2]+'\t'+str(pos6[0])+'\t'+str(pos4[1])+'\n')
	#Now it's time to detect inversion of the first or the second fragments of reciprocal translocation
	if TYPE == '0':
		for n in dic_CHR_ff:
			for j in dic_CHR_ff[n]:
				pos1 = [j[0], j[1]]
				pos3 = [j[3], j[4]]
				pos2 = []
				pos4 = []
				if j[2] in dic_CHR_fr_rev:
					for k in dic_CHR_fr_rev[j[2]]:
						if pos3[0] <= k[1] + 2*INSERT and  pos3[0] + vasouille >= k[1] and n == k[2]:
							pos4 = [k[0], k[1]]
							pos2 = [k[3], k[4]]
							if n in dic_CHR_rr:
								for x in dic_CHR_rr[n]:
									if pos2[0] + vasouille >= x[1] and pos2[0] <= x[1] + 2*INSERT and j[2] == x[2]:
										pos5 = [x[0], x[1]]
										pos6 = [x[3], x[4]]
										if j[2] in dic_CHR_rf_rev:
											for y in dic_CHR_rf_rev[j[2]]:
												if y[0] + vasouille >= pos6[1] and y[0] <= pos6[1] + 2*INSERT and y[4] <= pos1[0] + vasouille and y[4] + 2*INSERT >= pos1[0] and n == y[2]:
													pos8 = [y[0], y[1]]
													pos7 = [y[3], y[4]]
													if pos2[0] + vasouille >= pos1[1]:#first fragment inversed
														outfile.write('reciprocal_translocation\tregion1_inv:\t'+n+'\t'+str(pos1[0])+'\t'+str(pos5[1])+'\tregion2:\t'+j[2]+'\t'+str(pos8[0])+'\t'+str(pos4[1])+'\n')
													elif pos2[1] <= pos1[0] + vasouille:#second fragment inversed
														outfile.write('reciprocal_translocation\tregion1:\t'+n+'\t'+str(pos2[0])+'\t'+str(pos7[1])+'\tregion2_inv:\t'+j[2]+'\t'+str(pos3[0])+'\t'+str(pos6[1])+'\n')
	outfile.close()

# def couv_med_reg(DEB, FIN, CHR, COUV):
	# file = open(COUV)
	# LISTE = []
	# for line in file:
		# data = line.split()
		# if data:
			# if data[0] == CHR:
				# if int(data[1]) >= DEB and int(data[1]) <= FIN:
					# LISTE.append(int(data[2]))
				# elif int(data[1]) > FIN:
					# break
	# if len(LISTE) == 0:
		# return 0
	# else:
		# return mediane(LISTE)
		
def couv_med_reg(deb, fin, chr):
	"""
		:param chr: The chromosome name
		:type chr: str
		:param deb: the start position
		:type deb: int
		:param end: The end position
		:type end: int
	"""
	subListNoGap = []
	for site in covDic[chr][deb+1:fin+2]:
		if site:
			subListNoGap.append(site)
	if len(subListNoGap) == 0:
		return 0
	else:
		return mediane(subListNoGap)

# def couv_med_reg_del(DEB, FIN, CHR, COUV):
	# file = open(COUV)
	# LISTE = []
	# for line in file:
		# data = line.split()
		# if data:
			# if data[0] == CHR:
				# if int(data[1]) >= DEB and int(data[1]) <= FIN:
					# LISTE.append(int(data[2]))
				# elif int(data[1]) > FIN:
					# break
	# if len(LISTE) <= 0.5*(FIN-DEB):
		# return 0
	# else:
		# return mediane(LISTE)

def couv_med_reg_del(deb, fin, chr):
	"""
		:param chr: The chromosome name
		:type chr: str
		:param deb: the start position
		:type deb: int
		:param end: The end position
		:type end: int
	"""
	subListNoGap = []
	for site in covDic[chr][deb+1:fin+2]:
		if site:
			subListNoGap.append(site)
	if len(subListNoGap) <= 0.5*(fin-deb):
		return 0
	else:
		return mediane(subListNoGap)
	

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

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script Try to identify signature for structural variation")
	# Wrapper options.
	parser.add_option( '', '--frf', dest='frf', default='not_filled', help='The discordant fr or rf file depending on expected orientation')
	parser.add_option( '', '--ff', dest='ff', default='not_filled', help='The discordant ff file')
	parser.add_option( '', '--rr', dest='rr', default='not_filled', help='The discordant rr file')
	parser.add_option( '', '--ins', dest='ins', default='not_filled', help='The discordant ins file')
	parser.add_option( '', '--delet', dest='delet', default='not_filled', help='The discordant del file')
	parser.add_option( '', '--chr_rr', dest='chr_rr', default='not_filled', help='The discordant chr_rr file')
	parser.add_option( '', '--chr_fr', dest='chr_fr', default='not_filled', help='The discordant chr_fr file')
	parser.add_option( '', '--chr_rf', dest='chr_rf', default='not_filled', help='The discordant chr_rf file')
	parser.add_option( '', '--chr_ff', dest='chr_ff', default='not_filled', help='The discordant chr_ff file')
	parser.add_option( '', '--chr', dest='chr', default='not_filled', help='The tabulated file containing in col 1 : chromosome name, col 2: chromosome length. A line for each chromosomes')
	parser.add_option( '', '--covf', dest='covf', default='not_filled', help='The coverage file calculated in "calc_stat"')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--insert', dest='insert', default=5000, help='The expected insert size, [default: %default]')
	parser.add_option( '', '--exp_cov', dest='exp_cov', default='not_filled', help='The expected coverage (float)')
	parser.add_option( '', '--ploid', dest='ploid', default=0.33, help='Multiplicator for coverage variation detection in SV identification (ex : If homozygous duplication expected in diploid: expected = coverage + coverage*1, if heterozygous duplication expected in diploid => expected = coverage + coverage*0.5). Choose a value lower than the expected one')
	parser.add_option( '', '--type', dest='type', default='not_filled', help='The type of SV searched: q')
	parser.add_option( '', '--out', dest='out', default='not_filled', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	
	if options.frf == 'not_filled':
		mot = 'Please provide an argument for --frf'
		sys.exit(mot)
	if options.ff == 'not_filled':
		mot = 'Please provide an argument for --ff'
		sys.exit(mot)
	if options.rr == 'not_filled':
		mot = 'Please provide an argument for --rr'
		sys.exit(mot)
	if options.rr == 'not_filled':
		mot = 'Please provide an argument for --rr'
		sys.exit(mot)
	if options.chr_rr == 'not_filled':
		mot = 'Please provide an argument for --chr_rr'
		sys.exit(mot)
	if options.chr_rf == 'not_filled':
		mot = 'Please provide an argument for --chr_rf'
		sys.exit(mot)
	if options.chr_fr == 'not_filled':
		mot = 'Please provide an argument for --chr_fr'
		sys.exit(mot)
	if options.chr_ff == 'not_filled':
		mot = 'Please provide an argument for --chr_ff'
		sys.exit(mot)
	if options.ins == 'not_filled':
		mot = 'Please provide an argument for --ins'
		sys.exit(mot)
	if options.chr_rr == 'not_filled':
		mot = 'Please provide an argument for --chr_rr'
		sys.exit(mot)
	if options.type == 'not_filled':
		mot = 'Please provide an argument for --type'
		sys.exit(mot)
	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		sys.exit(mot)
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		covDic = {}
		if options.orient == 'rf':
			# list2rm = couv2couv(config.get('General','chr'), config.get('Calc_coverage','out'), tmp)
			couv2couv(config.get('General','chr'), config.get('Calc_coverage','out'), covDic)
			indent_discord(options.ff, options.frf, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, config.getfloat('Calc_coverage', 'median_insert'), options.out, config.getfloat('Calc_coverage', 'median_coverage'), config.getfloat('General','ploid'), options.type)
		elif options.orient == 'fr':
			# list2rm = couv2couv(config.get('General','chr'), config.get('Calc_coverage','out'), tmp)
			couv2couv(config.get('General','chr'), config.get('Calc_coverage','out'), covDic)
			indent_discord(options.rr, options.frf, options.ff, options.ins, options.delet, options.chr_ff, options.chr_rf, options.chr_fr, options.chr_rr, config.getfloat('Calc_coverage', 'median_insert'), options.out, config.getfloat('Calc_coverage', 'median_coverage'), config.getfloat('General','ploid'), options.type)
		else:
			mot = 'Unrecognized orientation: '+options.orient
			sys.exit(mot)
	else:
		if options.exp_cov == 'not_filled':
			mot = 'Please provide an argument for --exp_cov'
			sys.exit(mot)
		if options.covf == 'not_filled':
			mot = 'Please provide an argument for --covf'
			sys.exit(mot)
		tmp = (tempfile.NamedTemporaryFile().name)+'.cov'
		if options.orient == 'rf':
			list2rm = couv2couv(options.chr, options.covf, tmp)
			indent_discord(options.ff, options.frf, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, float(options.insert), options.out, tmp, float(options.exp_cov), float(options.ploid), options.type)
		elif options.orient == 'fr':
			list2rm = couv2couv(options.chr, options.covf, tmp)
			indent_discord(options.rr, options.frf, options.ff, options.ins, options.delet, options.chr_ff, options.chr_rf, options.chr_fr, options.chr_rr, float(options.insert), options.out, tmp, float(options.exp_cov), float(options.ploid), options.type)
		else:
			mot = 'Unrecognized orientation: '+options.orient
			sys.exit(mot)
	for n in list2rm:
		os.remove(n)

if __name__ == "__main__": __main__()
import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, random, multiprocessing

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


# def hold_job2(ID_liste):
	# time.sleep(10)
	# qs=os.popen("qstat")
	# nb_run = 0
	# for n in qs:
		# l = n.split()
		# if len(l) > 2:
			# if l[0] in ID_liste:
				# nb_run = nb_run + 1
	# while nb_run > 0:
		# time.sleep(10)
		# qs=os.popen("qstat")
		# nb_run = 0
		# for n in qs:
			# l = n.split()
			# if len(l) > 2:
				# if l[0] in ID_liste:
					# nb_run = nb_run + 1

# def job_ID2list(MOT):
	# for n in MOT:
		# return(n.split()[0])

# def cherche_error(FILE):
	# file = open(FILE)
	# trouve = 0
	# mot = ''
	# for line in file:
		# data = line.split()
		# if data:
			# trouve = 1
			# mot = mot + line
	# print mot


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n\n"
	"This script try to identify signature for structural variation (this is a wrapper for idend_sv.py)")
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
	parser.add_option( '', '--covf', dest='covf', default='not_filled', help='The coverage file calculated in "5_calc_stat"')
	parser.add_option( '', '--orient', dest='orient', default='rf', help='The expected orientation: rf or fr, [default: %default]')
	parser.add_option( '', '--insert', dest='insert', default=5000, help='The expected insert size, [default: %default]')
	parser.add_option( '', '--exp_cov', dest='exp_cov', default='not_filled', help='The expected coverage (float)')
	parser.add_option( '', '--ploid', dest='ploid', default=0.33, help='Multiplicator for coverage variation detection in SV identification (ex : If homozygous duplication expected in diploid: expected = coverage + coverage*1, if heterozygous duplication expected in diploid: expected = coverage + coverage*0.5). Choose a value lower than the expected one')
	parser.add_option( '', '--thread', dest='thread', default='1', help='The thread number used for search (integer), [default: %default]')
	parser.add_option( '', '--out', dest='out', default='SV_detected.tab', help='Output file')
	parser.add_option( '', '--config', dest='config', default=None)
	(options, args) = parser.parse_args()
	
	
	
	proc = int(options.thread)
	
	ScriptPath = os.path.dirname(sys.argv[0])
	
	loca_programs = ConfigParser.RawConfigParser()
	loca_programs.read(ScriptPath+'/loca_programs.conf')
	
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
	if options.out == 'not_filled':
		mot = 'Please provide an argument for --out'
		sys.exit(mot)
	i = 0
	liste_tmp = []
	liste_id = []
	# liste_job = []
	while i < 10:
		temp = options.out+'_'+str(i)
		liste_tmp.append(temp)
		# print temp
		if options.config:
			# qs=os.popen('qsub -q bioinfo.q -terse -b yes -V -N IDENT_SV "python %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --config %s --type %s"' % (ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, options.config, str(i)))
			# liste_job.append(job_ID2list(qs))
			liste_id.append("%s %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --config %s --type %s" % (loca_programs.get('Programs','python'), ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, options.config, str(i)))
		else:
			# qs=os.popen('qsub -q bioinfo.q -terse -b yes -V -N IDENT_SV "python %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --type %s"' % (ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, str(i)))
			# liste_job.append(job_ID2list(qs))
			liste_id.append("%s %s/ident_SV.py --frf %s --ff %s --rr %s --ins %s --delet %s --chr_rr %s --chr_fr %s --chr_rf %s --chr_ff %s --chr %s --covf %s --orient %s --insert %s --exp_cov %s --ploid %s --out %s --type %s" % (loca_programs.get('Programs','python'), ScriptPath, options.frf, options.ff, options.rr, options.ins, options.delet, options.chr_rr, options.chr_fr, options.chr_rf, options.chr_ff, options.chr, options.covf, options.orient, options.insert, options.exp_cov, options.ploid, temp, str(i)))
		i += 1
	
	liste_process = []
	for n in liste_id:
		t = multiprocessing.Process(target=run_job, args=(n, 'Bug lauching indent_SV.py',))
		liste_process.append(t)
		if len(liste_process) == proc:
			# Starts threads
			for process in liste_process:
				process.start()
			# This blocks the calling thread until the thread whose join() method is called is terminated.
			for process in liste_process:
				process.join()
			#the processes are done
			liste_process = []
	if liste_process:
		# Starts threads
		for process in liste_process:
			process.start()
		# This blocks the calling thread until the thread whose join() method is called is terminated.
		for process in liste_process:
			process.join()
		#the processes are done
		liste_process = []
	
	if options.config:
		config = ConfigParser.RawConfigParser()
		config.read(options.config)
		config.set('Ident_discord','frf', options.frf)
		config.set('Ident_discord','ff', options.ff)
		config.set('Ident_discord','rr', options.rr)
		config.set('Ident_discord','ins', options.ins)
		config.set('Ident_discord','delet', options.delet)
		config.set('Ident_discord','chr_rr', options.chr_rr)
		config.set('Ident_discord','chr_fr', options.chr_fr)
		config.set('Ident_discord','chr_rf', options.chr_rf)
		config.set('Ident_discord','chr_ff', options.chr_ff)
		with open(options.config, 'wb') as configfile:
			config.write(configfile)
	# for n in liste_job:
		# cherche_error('IDENT_SV.o'+n)
		# os.system('rm IDENT_SV.o'+n)
	mot = 'cat '
	for n in liste_tmp:
		mot = mot + n + ' '
	mot = mot + '> ' + options.out
	os.system(mot)
	for n in liste_tmp:
		os.remove(n)
		
if __name__ == "__main__": __main__()
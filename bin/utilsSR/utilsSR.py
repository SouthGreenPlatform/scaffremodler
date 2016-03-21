import optparse, os, shutil, subprocess, sys, tempfile, fileinput, ConfigParser, operator, time, math, glob, datetime
from inspect import currentframe, getframeinfo

def stop_err( msg ):
	sys.stderr.write( "%s\n" % msg )
	raise SystemExit()

def run_job (frameinfo, cmd_line, ERROR, silent = False):

	if not silent:
		print (cmd_line)
	try:
		tmp = tempfile.NamedTemporaryFile().name
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
		stop_err( 'Line : '+str(frameinfo.lineno)+' - '+ERROR + str( e ) )

def run_job_silent (frameinfo, cmd_line, ERROR):
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
		stop_err( 'Line : '+str(frameinfo.lineno)+'\n'+cmd_line+'\n'+ERROR + str( e ) )

def extractSamFromPosition(LOCA_PROGRAMS, SAM, TYPE, CHR, START, END, OUT):
	"""
		Provide a sam file or bam file from coordinates

		This function create a sam or bam file, according to the input format, from a first sam or bam file and coordinates like chr01:10000-20000.

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam or bam file
		:type SAM: str
		:param TYPE: The format of the input file. If the input file is a bam file, he must be indexed.
		:type TYPE: str ("sam" | "bam")
		:param CHR: File containing col1: chromosome name and col2: chromsome size
		:type CHR: str
		:param START: The start position
		:type START: int
		:param END: The end position
		:type END: int
		:param OUT: The name of the output file
        :return: void
	"""

	if TYPE == 'bam':
		bam2subbam = '%s view -bh %s %s:%s-%s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, CHR, START, END, OUT)
		run_job(getframeinfo(currentframe()), bam2subbam, 'Error in bam2subbam:\n', True)
	elif TYPE == 'sam':
		TMP = tempfile.NamedTemporaryFile().name.split('/')[-1]
		outbedTempName = TMP+'_outbed_tmp'
		outbed = open(outbedTempName, 'w')
		outbed.write(CHR+'\t'+START+'\t'+END)
		outbed.close()

		bam2subbam = '%s view -Sh -L %s -o %s %s' % (LOCA_PROGRAMS.get('Programs','samtools'), outbedTempName, SAM, OUT)
		run_job(getframeinfo(currentframe()), bam2subbam, 'Error in bam2subbam:\n', True)
		os.remove(outbedTempName)
	else:
		raise ValueError('File type is invalid.')


def indexBamFile(LOCA_PROGRAMS, BAM):
	"""
		Index a bam file.

		The name of the index file is : input_name.bai


		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param BAM: The input bam file.
		:return: void
	"""
	indexBamFile = '%s index %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM)
	run_job(getframeinfo(currentframe()), indexBamFile, 'Error in index_bam_file:\n', True)


def sam2bam(LOCA_PROGRAMS, SAM, OUT):
	"""
		Convert a sam file to bam format

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam file.
		:type SAM: str
		:param OUT: The output file name
		:type OUT: str
		:return: void
	"""
	sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	run_job(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam:\n', True)


def bam2sam(LOCA_PROGRAMS, BAM, OUT):
	"""
		Convert a bam file to sam format

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param BAM: The input bam file.
		:type BAM: str
		:param OUT: The output file name
		:type OUT: str
		:return: void
	"""
	bam2sam = '%s view -h %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), BAM, OUT)
	run_job(getframeinfo(currentframe()), bam2sam, 'Error in bam conversion to sam:\n', True)


def extractSamFromReads(locaPrograms, bam, type, READS, OUT):
	"""
		Provide a sam file from a bam or sam file containing only alignments corresponding to a reads list

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param bam: The input bam | sam file.
		:type bam: str
		:param type: The format of the input file ("sam" | "bam")
		:type type: str
		:param READS: A one column file containing the name of reads to extractSamFromPosition
		:type READS: str
		:param OUT: The output file name
		:type OUT: str
		:return: void
	"""

	if type == 'sam':
		tmp = tempfile.NamedTemporaryFile().name.split('/')[-1]
		sam2bam(locaPrograms, bam, tmp)
		FileToExtract = tmp
	elif type == 'bam':
		FileToExtract = bam
	else:
		raise ValueError("Format inconnu")

	#extract the header of the bam file
	recupHeader = '%s view -H %s -o %s' % (locaPrograms.get('Programs','samtools'), FileToExtract, OUT)
	run_job(getframeinfo(currentframe()), recupHeader, 'Error in extractSamFromReads :\n', True)

	#concatenate the alignments specified to the header
	extractReads = '%s -R %s %s >> %s' % (locaPrograms.get('Programs','bamgrepreads'), READS, FileToExtract, OUT)
	run_job(getframeinfo(currentframe()), extractReads, 'Error in extractSamFromReads:\n', True)

	#remove the temporary converted bam file
	if type == 'sam':
		os.remove(FileToExtract)


def calcul_cov(LOCA_PROGRAMS, SAM, TYPE, OUT):
	"""
		Calculate the coverage of a sam or bam file, site by site

		:param LOCA_PROGRAMS: From the Configparser module. Contains the path of each programs
		:param SAM: The input sam or bam file.
		:type SAM: str
		:param TYPE: The format of the input file
		:type TYPE: str ("sam" | "bam")
		:param OUT: The name of the output file
		:type OUT: str
        :return: void
	"""

	if TYPE == 'sam':
		#Convert the sam file to the bam format
		sam2bam = '%s view -bSh %s -o %s' % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, SAM+'_BAM.bam')
		run_job(getframeinfo(currentframe()), sam2bam, 'Error in sam2bam (calcul_cov):\n', True)
		SAM = SAM+'_BAM.bam'
	elif TYPE != 'bam':
		mot = TYPE+' argument passed in --type is not recognized'
		raise ValueError(mot)

	cal_cov = ('%s depth %s > %s') % (LOCA_PROGRAMS.get('Programs','samtools'), SAM, OUT)
	run_job(getframeinfo(currentframe()), cal_cov, 'Error in calculating coverage:\n', True)

	#remove the intermediate bam file.
	if TYPE == 'sam':
		os.remove(SAM)


def mediane(L):
	"""
		Give the median value of a list of integers

		:param L: A list of integers
		:type L: list
		:return: The median value
        :rtype: int
	"""

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

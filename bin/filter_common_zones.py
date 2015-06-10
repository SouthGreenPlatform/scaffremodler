
#
#  Copyright 2014 CIRAD
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

# -*- coding: utf-8 -*-
import ConfigParser
import datetime
import math
import multiprocessing
from multiprocessing import Queue
import optparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import utilsSR.utilsSR as utils
from inspect import currentframe, getframeinfo
from operator import itemgetter


def checkFileExist(fname):
	"""
		Check the existence of a file.

		:param fname: the path to the file to checkFileExist
		:type fname: str
		:return: A boolean with true if the file existe, or false if the file doesn't exist.
		:rtype: bool
	"""

	try:
		f = open(fname,'r')
		f.close()
		return 1
	except:
		return 0


def parseConfFile(confFile, locaPrograms):
	"""
		Validate the configuration file, with no duplicate sections or options in a section.

		:param confFile: Path to the configuration file
		:type confFile: str
		:param locaPrograms: Path to the file containing all the tools paths.
		:type locaPrograms: str
		:return: A list of 2 elements. First : true for a validate conf, or false otherwise.
				 Second : The error message if there is an error in the conf file. Empty element ortherwise.
		:rtype: list
	"""

	if checkFileExist(confFile):

		listJobs = []
		discType = ""
		job = {}
		f = open(confFile, 'r')
		for line in f:
			line = line.strip()
			if line.strip():
				if line[0] == '[':
					if job:
						listJobs.append([job, locaPrograms])
					job = {}
					discType = line[1:-1]
					if not discType in job:
						job[discType] = {}
					else:
						raise ValueError("Error in the configuration file\nDuplicate section \""+discType+"\" found.")
				else:
					elmnts = re.split("[=:]", line)
					if not elmnts[0] in job[discType]:
						genomeID = elmnts[0].strip()
						job[discType][genomeID] = {}
						elmts = elmnts[1].strip().split()

						job[discType][genomeID] = {'scoreFile': elmts[0], \
												   'bamFile': elmts[1], \
												   'minReadsNumber': int(elmts[2]), \
												   'libraryLength': int(elmts[3]), \
												   'margin': int(elmts[4])
												  }
					else:
						raise ValueError("Error in the configuration file\nDuplicate option \""+elmnts[0]+"\" found in the \""+discType+"\" section.")
		if job:
			listJobs.append([job, locaPrograms])
	else:
		raise IOError("Configuration file : " +confFile+" is not found.")

	return listJobs


def defineFlag(tag):
	"""
		Create a flag from a tag.

		For tag "NEW_PASSED", create a flag as "1"
		For tag "NOT_PASSED", create a flag as "2"
		For tag "PASSED", create a flag as "3"

		:param tag: the tag...
		:type tag: str
		:return: a string corresponding to the flag
		:rtype: str
	"""
	return {
        'NEW_PASSED': "1",
        'NOT_PASSED': "2",
		'PASSED': "3"
        }[tag]


def zonesOverlap(zonesOrg, zonesCompList):
	"""
		Calculate the overlap of one couple of zone with all the zones of one accession.

		:param zonesOrg: The mate zone to test.
		:type zonesOrg: list
		:param zonesCompList: The zones of one accession
		:type zonesCompList: list
		:return: A list with three boolean. first : indicate if the zoneOrg is overlapping one zone of zonesCompList; The second indicate if the overlapping zone is tagged as "PASSED" or not. The third indicate if the overlapping zone is tagged as "NOT_PASSED" or not.
		:rtype: list
	"""
	i = 0
	j = len(zonesCompList)
	m = (i + j) // 2

	overlap = False
	zonePassed = False
	zoneNotPassed = False

	amontChr = zonesOrg[0]
	amontStart = int(zonesOrg[1])
	amontEnd = int(zonesOrg[2])

	avalChr = zonesOrg[5]
	AvalStart = int(zonesOrg[6])
	AvalEnd = int(zonesOrg[7])

	while i < j and not overlap and not zonePassed:

		if zonesCompList[m][0] == amontChr:

			k = m
			while k >= i and not zonePassed and zonesCompList[k][0] == amontChr and zonesCompList[k][5] == avalChr:

				if amontStart <= int(zonesCompList[k][2]) and amontEnd >= int(zonesCompList[k][1]) and AvalStart <= int(zonesCompList[k][7]) and AvalEnd >= int(zonesCompList[k][6]):
					overlap = True
					if zonesCompList[m][-1] == "PASSED":
						zonePassed = True
					elif zonesCompList[m][-1] == "NOT_PASSED":
						zoneNotPassed = True
				k -= 1

			k = m + 1
			while k < j and not zonePassed and zonesCompList[k][0] == amontChr and zonesCompList[k][5] == avalChr:

				if amontStart <= int(zonesCompList[k][2]) and amontEnd >= int(zonesCompList[k][1]) and AvalStart <= int(zonesCompList[k][7]) and AvalEnd >= int(zonesCompList[k][6]):
					overlap = True
					if zonesCompList[m][-1] == "PASSED":
						zonePassed = True
					elif zonesCompList[m][-1] == "NOT_PASSED":
						zoneNotPassed = True
				k += 1

			# no infinite loop
			if not overlap:
				break

		elif zonesCompList[m][0] < amontChr:
			i = m + 1
		else:
			j = m

		m = (i + j)//2

	return [overlap, zonePassed, zoneNotPassed]


def amont_aval(LocalProgram, sam, type, output):
	"""
		Create two sam files with separate alignments according to the position of the reads.

		Parse a sam or bam file witch contains pair end alignments and separate the reads according to their position.
		When the reads are aligned on two different chromosomes, the order of the chromosomes is given by the parameter CHR.

		:param LocalProgram: From the Configparser module. Contains the path of each programs
		:param sam: The input sam or bam file.
		:type sam: str
		:param type: The format of the input file. If the input file is a bam file, he must be indexed.
		:type type: str ("sam" | "bam")
		:param output: The name of the output file
		:type output: str
        :return: void

		.. seealso:: sort_sam()
        .. warnings:: This function only work in paired end reads.
	"""
	if type == 'bam':

		tmpOutput = tempfile.NamedTemporaryFile(delete=False, dir=os.getcwd())
		tmpOutput.close()
		utils.bam2sam(LocalProgram, sam, tmpOutput.name)
		total = sort_sam(tmpOutput.name, output)
		os.remove(tmpOutput.name)
		return total
	elif type == 'sam':
		total = sort_sam(sam, output)
		return total
	else:
		sys.exit('Unrecognised argument passed in --type option')


def sort_sam(sam, output):

	"""
		Create two sam files with separate alignments according to the position of the reads.

		Parse a sam or bam file witch contains pair end alignments and separate the reads according to their position.
		When the reads are aligned on two different chromosomes, the order of the chromosomes is given by the parameter CHR.

		:param sam: The input sam or bam file.
		:type sam: str
		:param output: The name of the output file
		:type output: str
        :return: The number of common reads in the two mate zones
		:rtype: int

		.. seealso:: amont_aval()
	"""

	#parsing sam file
	amont = open(output+'_amont_sam','w')
	aval = open(output+'_aval_sam','w')
	file = open(sam)
	total = 0
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == '@': # lines corresponding to the sam header's
				amont.write(line)
				aval.write(line)
			else:
				if data[6] == '=' or data[6] == data[2]:
					if int(data[3]) <= int(data[7]):
						total += 1
						amont.write(line)
					elif int(data[3]) > int(data[7]):
						total += 1
						aval.write(line)
				else:
					if data[2] < data[6]:
						amont.write(line)
					else:
						aval.write(line)
					total += 1

	amont.close()
	aval.close()
	return total/2


def concatAllZones(listZones, filter=''):
	"""
		Concatenate the zones from several lists

		:param listZones: A list containing several list of zones.
		:type listZones: list
		:param filter: A filter to merge only the zone where the tag in the 13th column corresponding to the filter.
		:type filter: str
		:return: A list with all the zones concatenated.
		:rtype: list
	"""

	allZones = []
	for genomeZones in listZones:
		for zone in listZones[genomeZones]:
			if filter and zone[13] != filter:
				continue
			else:
				allZones.append(zone)
	return allZones


def mergeZones(locaPrograms, listZones, filter=''):
	"""
		Merge overlapping zones.

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param listZones: A list containing the zones to merge.
		:type listZones: list
		:param filter: A filter to merge only the zone where the tag in the 13th column corresponding to the filter.
		:type filter: str
		:return: A list with 2 elements. The fist element is a list with all the zone non merged and/or filtered. THe second is a list with the zone merged.
		:rtype: list
	"""

	allZones = []  # containing all the zones to try to merge, based on the filter
	zonesOrg = []  # containing the no merged zones and the filtered zones
	zonesMerged = []  # containing the zones merged

	if filter:
		for zone in listZones:
			if zone[13] == filter:
				allZones.append(zone)
			else:
				zonesOrgs.append(zone)
	else:
		allZones = listZones

	d = sorted(allZones, key=itemgetter(1))
	allZones = sorted(d, key=itemgetter(0))

	i = 0
	zoneFound = False
	for zone in allZones:
		if zone:
			k = i - 1
			while k > 0:
				zoneCat = zone
				if allZones[k]:
					if allZones[k][0] == zoneCat[0] and int(allZones[k][1]) <= int(zoneCat[2]) and int(allZones[k][2]) >= int(zoneCat[1]):
						if allZones[k][5] == zoneCat[5] and int(allZones[k][6]) <= int(zoneCat[7]) and int(allZones[k][7]) >= int(zoneCat[6]):
							amontPosStart = min(int(allZones[k][1]), int(zoneCat[1]))
							amontPosEnd = max(int(allZones[k][2]), int(zoneCat[2]))

							avalPosStart = min(int(allZones[k][6]), int(zoneCat[6]))
							avalPosEnd = max(int(allZones[k][7]), int(zoneCat[7]))

							if len(allZones[k]) == 15:
								a = allZones[k][14]+zoneCat[14]
								zoneCat = [zoneCat[0], amontPosStart, amontPosEnd, '', '', zoneCat[5], avalPosStart, avalPosEnd, '', '', '', '', '', filter, a]
							else:
								zoneCat = [zoneCat[0], amontPosStart, amontPosEnd, '', '', zoneCat[5], avalPosStart, avalPosEnd, '', '', '', '', '', filter]

							allZones[k] = 0
							allZones[i] = zoneCat
							zoneFound = True
					else:
						break
				k -= 1

			k = i + 1
			while k < len(allZones):  # Try to find downstream overlapping zones
				zoneCat = zone
				if allZones[k]:
					if allZones[k][0] == zoneCat[0] and int(allZones[k][1]) <= int(zoneCat[2]) and int(allZones[k][2]) >= int(zoneCat[1]):
						if allZones[k][5] == zoneCat[5] and int(allZones[k][6]) <= int(zoneCat[7]) and int(allZones[k][7]) >= int(zoneCat[6]):
							amontPosStart = min(int(allZones[k][1]), int(zoneCat[1]))
							amontPosEnd = max(int(allZones[k][2]), int(zoneCat[2]))

							avalPosStart = min(int(allZones[k][6]), int(zoneCat[6]))
							avalPosEnd = max(int(allZones[k][7]), int(zoneCat[7]))

							if len(allZones[k]) == 15:
								a = allZones[k][14]+zoneCat[14]
								zoneCat = [zoneCat[0], amontPosStart, amontPosEnd, '', '', zoneCat[5], avalPosStart, avalPosEnd, '', '', '', '', '', filter, a]
							else:
								zoneCat = [zoneCat[0], amontPosStart, amontPosEnd, '', '', zoneCat[5], avalPosStart, avalPosEnd, '', '', '', '', '', filter]

							allZones[k] = 0
							allZones[i] = zoneCat
							zoneFound = True
					else:
						break
				k += 1

			if not zoneFound:
				allZones[i] = 0
				zonesOrg.append(zone)
		i += 1

	for zone in allZones:
		if zone:

			zonesMerged.append(zone)

	# Sort the list on chromosomes and first position of the upstream zone
	d = sorted(zonesOrg, key=itemgetter(1))
	zonesOrg = sorted(d, key=itemgetter(0))

	d = sorted(zonesMerged, key=itemgetter(1))
	zonesMerged = sorted(d, key=itemgetter(0))

	return [zonesOrg, zonesMerged]


def findZoneInBam(locaPrograms, bamAmont, bamAval, zones, nbReadsMax, margin):
	"""
		From mate zones and a bam file, find if the mate zones have commons reads in the bam file.

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param bamAmont: The bam file of upstream reads.
		:type bamAmont: str
		:param bamAval: The bam file of downstream reads.
		:type bamAval: str
		:param zones: The mate zone to test.
		:type zones: list
		:param nbReadsMax: The number of reads to validate the zones.
		:type nbReadsMax: int
		:param margin: The margin added to the border of the zone.
		:type margin: int
		:return: a list with the zone found from the bam files. If the zones are not validate with the nbReadsMax, return an empty list.
		:rtype: list
	"""
	zone = []  # to return
	utils.extractSamFromPosition(locaPrograms, bamAmont, 'bam', zones[0], str(int(zones[1])-int(margin)), str(int(zones[2])+int(margin)), 'tmp_amont_subsam')
	utils.extractSamFromPosition(locaPrograms, bamAval, 'bam', zones[3], str(int(zones[4])-int(margin)), str(int(zones[5])+int(margin)), 'tmp_aval_subsam')

	nbReadsCommons = intersectReads(locaPrograms, 'tmp_amont_subsam', 'tmp_aval_subsam')

	if nbReadsCommons > nbReadsMax:
		zoneAmont = calculBorder(locaPrograms, 'tmp_amont_subsam')
		zoneAval = calculBorder(locaPrograms, 'tmp_aval_subsam')

		zone = [zoneAmont[0],int(zoneAmont[1]),int(zoneAmont[2]),int(zoneAmont[2])-int(zoneAmont[1]),float(zoneAmont[3]), \
				zoneAval[0],zoneAval[1],zoneAval[2],int(zoneAval[2])-int(zoneAval[1]),float(zoneAval[3]),"-",nbReadsCommons,"-","NEW_PASSED"]

	return zone


def intersectReads (locaPrograms, samAmont, samAval):
	"""
		Count the number of common reads beetween two zones

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param samAmont: The sam file of upstream alignments
		:type samAmont: str
		:param samAval: The sam file of downstream alignments
		:type samAval: str
		:return: The number of common reads
		:rtype: int
	"""

	extractReads = locaPrograms.get('Programs','samtools')+' view '+samAmont
	p = subprocess.Popen(extractReads,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
	readsAmont, errors = p.communicate()

	nbReads = 0

	if readsAmont:
		extractReads = locaPrograms.get('Programs','samtools')+' view '+samAval
		p = subprocess.Popen(extractReads,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
		readsAval, errors = p.communicate()

		if readsAval:
			listReadsAmont = []
			listReadsAval = []

			listAlignmentsAmont = readsAmont.split('\n')
			for alignments in listAlignmentsAmont:
				if alignments.strip():
					elements = alignments.split('\t')
					listReadsAmont.append(elements[0])

			listAlignmentsAval = readsAval.split('\n')
			for alignments in listAlignmentsAval:
				if alignments.strip():
					elements = alignments.split('\t')
					listReadsAval.append(elements[0])

			# intersect of the two zones list of reads
			setListAval = set(listReadsAval)
			intersect = [val for val in listReadsAmont if val in setListAval]
			nbReads = len(intersect)

	return nbReads


def calculBorder(locaPrograms, sam):
	"""
		calculate the border of zone, from the coverage of the zone.

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param sam: The sam file containing the alignments to calculate the borders.
		:type sam: str
		:return: A list like [chromosome, start_position, end_position, median_coverage]
		:rtype: list
	"""

	utils.calcul_cov(locaPrograms, sam, 'sam', 'tmp_sub_cov')

	cov = open('tmp_sub_cov', 'r')
	chrName = ""
	posStart = ""
	posEnd = ""
	listCovPos = []
	for line in cov:
		if line.strip():
			elments = line.split()
			if chrName == "":
				chrName = elments[0]
				posStart = elments[1]
				posEnd = elments[1]
				listCovPos.append(int(elments[2]))
			elif elments[0] == chrName:
				chrName = elments[0]
				posEnd = elments[1]
				listCovPos.append(int(elments[2]))
			else:
				raise ValueError("Error in calculBorder() : More than one chromosome found in a zone.")
	covMedian = utils.mediane(listCovPos)
	os.remove('tmp_sub_cov')
	return [chrName, posStart, posEnd, covMedian]


def calcAllZonesFromBam(locaPrograms, genomesList, listZones):
	"""
		From a bam file and a list of zone, extract informations about the zones

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param genomesList: A dictionary containing the information about the genomes
		:type genomesList: dict
		:param listZones: A list of all zones to extract
		:type listZones: list
		:return: void
	"""

	# Create a file with all the zones from each genomes merged.
	zonesConcat = concatAllZones(listZones, "PASSED")

	mergedZones = mergeZones(locaPrograms, zonesConcat)

	allZonesMerged = mergedZones[0]+mergedZones[1]
	del mergedZones, zonesConcat

	newZones = []
	for genome in genomesList:
		if genomesList[genome]['scoreFile'] == '-':
			listZones[genome] = []
			for zone in allZonesMerged:
				if zone:
					rslt = findZoneInBam(locaPrograms,
										 genome+'_amont_bam',
										 genome+'_aval_bam',
										 [zone[0], zone[1], zone[2], zone[5], zone[6], zone[7]],
										 genomesList[genome]['minReadsNumber'],
										 genomesList[genome]['margin'])
					if rslt:
						newZones.append(rslt)

			# sort the list on the first and the second elements
			d = sorted(newZones, key=itemgetter(1))
			f = sorted(d, key=itemgetter(5))
			newZones = sorted(f, key=itemgetter(0))

			output = open(genome+'.score', 'w')
			for zone in newZones:
				output.write('\t'.join(map(str, zone))+'\n')
				listZones[genome].append(zone)
			output.close()
			genomesList[genome]['scoreFile'] = genome+'.score'


def lookForCommonZones(locaPrograms, genomesList, discType):
	"""
		Search for common zones in all genomes in pairs.

		:param locaPrograms: From the Configparser module. Contains the path of each programs
		:param genomesList: A dictionary containing the information about the genomes
		:type genomesList: dict
		:param discType: The name of the discordant type.
		:type discType: str
	"""

	listZones = {}  # A dict containing a list of each zones by genome
	listNewZones = {}  # A dict containing a list of each new zone by genome
	noScoreFile = False

	for genome in genomesList:

		scoreFile = genomesList[genome]['scoreFile']
		bamFile = genomesList[genome]['bamFile']
		minReadsNumber = genomesList[genome]['minReadsNumber']
		libraryLength = genomesList[genome]['libraryLength']

		listZones[genome] = []  # Stock the zones which were described in the score file
		listNewZones[genome] = []  # Stock the new zones, obtained from the bam file.

		#index the bam file
		utils.indexBamFile(locaPrograms, bamFile)

		# sort the bam file with upstream and downstream alignments
		amont_aval(locaPrograms, genomesList[genome]['bamFile'], 'bam', genome)
		utils.sam2bam(locaPrograms, genome+'_amont_sam', genome+'_amont_bam')
		utils.indexBamFile(locaPrograms, genome+'_amont_bam')
		utils.sam2bam(locaPrograms, genome+'_aval_sam', genome+'_aval_bam')
		utils.indexBamFile(locaPrograms, genome+'_aval_bam')

		# Load in memory the zones of each genomes
		if scoreFile == '-':
			listZones[genome] = []
			noScoreFile = True
		else:
			if checkFileExist(scoreFile):
				f = open(scoreFile)
				listZones[genome] = []
				for line in f:
					if line.strip() and line[0] != '#':
						elmts = line.split()
						listZones[genome].append([elmts[0],int(elmts[1]),int(elmts[2]),int(elmts[3]),float(elmts[4]),elmts[5],int(elmts[6]),int(elmts[7]),int(elmts[8]),float(elmts[9]),elmts[10],int(elmts[11]),float(elmts[12]),elmts[13]])
				f.close()

				# sort the list on the first and the second elements
				d = sorted(listZones[genome], key=itemgetter(1))
				f = sorted(d, key=itemgetter(5))
				listZones[genome] = sorted(f, key=itemgetter(0))
			else:
				raise IOError("The score file "+scoreFile+" is not found.")

	# search for zones directly in a bam file if no score file found for one genome
	if noScoreFile:
		calcAllZonesFromBam(locaPrograms, genomesList, listZones)

	nbZones= {}

	for genome in genomesList:
		nbZones[genome] = {}
		i = 0

		for mateZones in listZones[genome]:
			if mateZones:
				listZones[genome][i].extend([[ ]])  # define a new element to store all the comparisons with each accessions
				for genomeComp in genomesList:
					if genomeComp != genome:
						zonesFound = zonesOverlap(mateZones, listZones[genomeComp])
						if zonesFound[0]:
							if zonesFound[1]:
								flag = "3"
							else:
								if zonesFound[2]:
									flag = "2"
								else:
									flag = "1"
							listZones[genome][i][14].append([genomeComp, flag])
							if genomeComp in nbZones[genome]:
								nbZones[genome][genomeComp] += 1
							else:
								nbZones[genome][genomeComp] = 1
						else:  # try to find reads directly in the bam file.
							zoneCoord = [mateZones[0], mateZones[1], mateZones[2], mateZones[5], mateZones[6], mateZones[7]]
							rslt = findZoneInBam(locaPrograms,
												 genomeComp+'_amont_bam',
												 genomeComp+'_aval_bam',
												 zoneCoord, genomesList[genome]['minReadsNumber'],
												 genomesList[genome]['margin'])
							if rslt:
								listZones[genome][i][14].append([genomeComp, "1"])
								flag = defineFlag(listZones[genome][i][13])
								rslt.extend([[ ]])
								rslt[14].append([genome, flag])
								listNewZones[genomeComp].append(rslt)
							else:
								listZones[genome][i][14].append([genomeComp, "0"])

			i += 1

	# write all the results in outfiles
	listOutput = {}
	parentDirPath = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

	for genome in genomesList:

		# sort the zones to be available to merge
		d = sorted(listNewZones[genome], key=itemgetter(1))
		newZones = sorted(d, key=itemgetter(0))

		newZonesMerged = mergeZones(locaPrograms, newZones)
		newZonesComp = []

		for zone in newZonesMerged[1]:
			rslt = findZoneInBam(locaPrograms,
								 genome+'_amont_bam',
								 genome+'_aval_bam',
								 [zone[0], zone[1], zone[2], zone[5], zone[6], zone[7]],
								 genomesList[genome]['minReadsNumber'],
								 genomesList[genome]['margin'])
			if rslt:
				rslt.append(zone[14])
				newZonesComp.append(rslt)

		allZones = listZones[genome]+newZonesMerged[0]+newZonesComp
		d = sorted(allZones, key=itemgetter(1))
		allZones = sorted(d, key=itemgetter(0))
		listOutput[genome] = open(genome+'_'+discType+'.score', 'w')

		for zone in allZones:
			genComp = {}
			flag = ""
			for comp in zone[14]:
				if comp[0] not in genComp:
					genComp[comp[0]] = comp[1]
				else:
					if comp[1] > genComp[comp[0]]:
						genComp[comp[0]] = comp[1]

			for genomeID in genomesList:
				if not genomeID in genComp and genomeID != genome:
					genComp[genomeID] = "0"

			for comp in genComp:
				if comp != genome:
					if flag:
						flag += "|"
					flag += comp + ':' + genComp[comp]

			listOutput[genome].write( \
									 zone[0]+'\t'+str(zone[1])+'\t'+str(zone[2])+'\t'+str(zone[3])+ '\t'+str(zone[4])+'\t'+ \
									 zone[5]+'\t'+str(zone[6])+'\t'+str(zone[7])+'\t'+str(zone[8])+'\t'+str(zone[9])+'\t'+ \
									 flag+'\t'+str(zone[11])+'\t'+str(zone[12])+'\t'+zone[13]+'\n')
		listOutput[genome].close()

		# move the score file to the parent directory
		os.rename(genome+'_'+discType+'.score', parentDirPath+'/'+genome+'_'+discType+'.score')

	# print ("\n##########\nResults for "+discType)
	# print(nbZones)
	# print("##########\n")
	return 0

def main(job):

	discType =  job[0].keys()[0]
	locaPrograms = job[1]

	# Create a dir based on the discordant type
	WORKING_DIR = tempfile.mkdtemp(prefix=discType+'_', dir=os.getcwd()).split('/')[-1]+'/'
	os.chdir(WORKING_DIR)
	try:
		rslt = lookForCommonZones(locaPrograms, job[0][discType], discType)
	except Exception as e:
		print e
		rslt = 1
	finally:
		os.chdir('../')
		shutil.rmtree(WORKING_DIR)
		return rslt


def worker(listJobs, out_q, locaPrograms):
	"""

	"""
	outdict = {}
	for job in listJobs:
		for discType in job:
			try:
				outdict[discType] = main(locaPrograms, job[discType], discType)
			except Exception as e:
				outdict[discType] = 0
				sys.stderr.write(format(e))
				pass
	out_q.put(outdict)

def __main__():

	t_start = datetime.datetime.now()

	# Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Paco Derouault : paco.derouault@gmail.com")
	parser.add_option( '-c', '--conf', dest='conf', default='', help='Config file')
	parser.add_option( '-p', '--nump', dest='nump', default='', help='Number of processor to use.')

	(options, args) = parser.parse_args()
	pathname = os.path.dirname(sys.argv[0])
	locaPrograms = ConfigParser.RawConfigParser()
	locaPrograms.read(pathname+'/loca_programs.conf')

	if not options.conf:
		mot = 'Please provide an argument for --gen'
		sys.exit(mot)

	if options.nump:
		nbProcs = int(options.nump)
	else:
		nbProcs = 1

	if nbProcs > multiprocessing.cpu_count():
		sys.exit("Processors number too high.\nYou have only "+str(multiprocessing.cpu_count())+" processor(s) available on this computer.")

	listJobs = parseConfFile(options.conf, locaPrograms)

	pool = multiprocessing.Pool(processes=nbProcs)
	results = pool.map(main, listJobs)

	# print(results)

	for i, job in enumerate(results):
		if job:
			print("Sorry the Discordant Type : "+listJobs[i][0].keys()[0]+" hasn't been completed.")


	# out_q = Queue()
	# chunksize = int(math.ceil(len(listJobs) / float(nbProcs)))  # Number of jobs per processor
	# procs = []


	# for i in range(nbProcs):
		# p = multiprocessing.Process(
				# target=worker,
				# args=(listJobs[chunksize * i:chunksize * (i + 1)],
					  # out_q,
					  # locaPrograms))
		# procs.append(p)
		# p.start()

	## Stock the results in a dict
	# resultdict = {}
	# for i in range(chunksize):
		# resultdict.update(out_q.get())

	# for p in procs:
		# p.join()

	# for job in resultdict:
		# if not resultdict[job]:
			# print("The job '"+job+"' could not be completed due to an error. Please read the error log file for more details.")

	print ("Total time : "+str(datetime.datetime.now()-t_start))

if __name__ == "__main__": __main__()

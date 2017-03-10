#template.vcf is needed to make a merged file
import operator
import vcf
import argparse
import sys

def extractInfoDelly(record):
	new_record = record
	idDelly = record.ID + "-DELLY"
	new_record.ID = idDelly
	new_record.REF = record.REF
	pos = record.POS
	end = record.INFO["END"]
	new_record.INFO["END"] = end
	svtype = record.INFO["SVTYPE"]

	if record.INFO["SVTYPE"] == "DEL":
		length = (end - pos) + 1
	elif record.INFO["SVTYPE"] == "DUP":
		length = (end - pos) + 1
	elif record.INFO["SVTYPE"] == "INV":
		length = (end - pos) + 1
	elif record.INFO["SVTYPE"] == "INS":
		length = record.INFO["INSLEN"]
	elif record.INFO["SVTYPE"] == "TRA":
		length = 0
	else:
		length = 0
		print "NO LENGTH", new_record.CHROM, pos, record.INFO["SVTYPE"]
	new_record.INFO["SVLEN"] = length

	return new_record


def extractInfoManta(record):
	new_record = record
	idManta = record.ID + "-MANTA"
	new_record.ID = idManta

	if "END" in record.INFO:
		end = record.INFO["END"]
	elif record.INFO["SVTYPE"] == "BND": #because BNDs do not have a value for END and contain mostly translocations with a length of 0, END is set equal to POS.
		end = record.POS
	else:
		print "there is a value missing for record.Info[\"END\"]"
		end = 0
		print record.CHROM, record.POS, record.INFO["SVTYPE"]
	new_record.INFO["END"] = end

	length = 0
	if "SVLEN" in record.INFO:
		length = record.INFO["SVLEN"]
		if type(length) != int:
			length = length[0]
	elif record.INFO["SVTYPE"] == "BND":
		length = 0
	elif record.INFO["SVTYPE"] == "INS":
		length = 0
		print record.INFO["SVTYPE"], record.CHROM, record.POS, "-length of insertion unknown, set to zero"
	new_record.INFO["SVLEN"] = length

	return new_record

def getInfoInList(record):
	new_record = record
	chrom = new_record.CHROM
	pos = new_record.POS
	end = new_record.INFO["END"]
	alt = new_record.ALT

	if "CIPOS" in new_record.INFO:
		cipos = new_record.INFO["CIPOS"]
		ciposmin = pos + cipos[0]
		ciposmax = pos + cipos[1]
	else:
		cipos = [-50, 50]
		ciposmin = pos + cipos[0]
		ciposmax = pos + cipos[1]

	if "CIEND" in new_record.INFO:
		ciend= new_record.INFO["CIEND"]
		ciendmin = end + ciend[0]
		ciendmax = end + ciend[1]
	else:
		ciend = [-50, 50]
		ciendmin = end + ciend[0]
		ciendmax = end + ciend[1]

	ciposint = abs(cipos[0]) + cipos[1] #the length of the CI
	ID = new_record.ID
	if "MANTA" in ID:
		tool = "MANTA"
	if "DELLY" in ID:
		tool = "DELLY"
	info = new_record.INFO
	form = new_record.FORMAT

	dict_record = {"ALT": alt, "ciposint": ciposint, "tool": tool, "CHROM" : chrom, "POS": pos, "ciposmin": ciposmin, "ciposmax": ciposmax, "ciendmin": ciendmin, "ciendmax": ciendmax, "ID" : ID, "INFO": info, "record": new_record}
	return dict_record

def startNewListForComparisonSVs(list_for_comparisonSVs, currentLine):
	comparingSVs_list = list_for_comparisonSVs
	if comparingSVs_list:
		comparingSVs_list = []
	comparingSVs_list.append(currentLine)

	return comparingSVs_list

def compareFilterSVs(list_for_comparisonSVs): #is called when list_for_comparisonSVs is about to be refreshed
	similar_SVs = list_for_comparisonSVs
	similar_SVs.sort(key=operator.itemgetter("ciposint"))
	similar_SVs.sort(key=operator.itemgetter("tool"), reverse = True)
	sv_to_print = []

	if len(similar_SVs) == 1:
		record_sv_to_print = similar_SVs[0]["record"]
		record_sv_to_print.INFO["CSA"] = 1
		sv_to_print.append(record_sv_to_print)
	else:
		if similar_SVs[0]["tool"] == "MANTA":
			record_sv_to_print = similar_SVs[0]["record"]
			for similar_SV in similar_SVs:
				if "DELLY" in similar_SV["tool"]:
					record_sv_to_print.INFO["CSA"] = 2
					record_sv_to_print.INFO["INFODELLY"] = similar_SV["POS"], #similar_SV["INFO"]["END"]
				else:
					record_sv_to_print.INFO["CSA"] = 1
			sv_to_print.append(record_sv_to_print)

		else: #Only delly
			record_sv_to_print = similar_SVs[0]["record"]
			record_sv_to_print.INFO["CSA"] = 1
			sv_to_print.append(record_sv_to_print)

	return sv_to_print

def filterBndAndWriteSVs(all_svs_to_write, vcf_writer):
	for record in all_svs_to_write:
		if record.INFO["SVTYPE"] == "BND":
			record_alt = record.ALT[0]
			string_record_alt = str(record_alt)
			if "[" in string_record_alt:
				split_string_record_alt = string_record_alt.split("[")
			if "]" in string_record_alt:
				split_string_record_alt = string_record_alt.split("]")
			splitted_alt = split_string_record_alt[1].split(":")
			chrom2 = splitted_alt[0]
			pos2 = splitted_alt[1]
			i_pos2 = int(pos2)
		else:
			chrom2 = 0; i_chrom2 = 0 # if no bnd, chrom should not match

		for line in all_svs_to_write:
			if (chrom2 == line.CHROM):
				if (i_pos2 == line.POS):
					all_svs_to_write.remove(line)

	for record in all_svs_to_write:
		vcf_writer.write_record(record)

def conditionsInsForComparison(currentLine, previousLine):
	#for insertions there is checked on length instead of END of SV
	svLenCurrentLine = currentLine["INFO"]["SVLEN"]
	svLenPreviousLine = previousLine["INFO"]["SVLEN"]
	cilength = [-20,20]
	svtypeCur = currentLine["INFO"]["SVTYPE"]
	svtypePrev = previousLine["INFO"]["SVTYPE"]
	ciposminCurrentLine = svLenCurrentLine + cilength[0]
	ciposmaxCurrentLine = svLenCurrentLine + cilength[1]
	ciposminPreviousLine = svLenPreviousLine + cilength[0]
	ciposmaxPreviousLine = svLenPreviousLine + cilength[1]

	if (((svtypeCur == "INS") and (svtypePrev == "INS")) and #for insertions
	(((ciposmaxCurrentLine >= ciposminPreviousLine) and (ciposminCurrentLine <= ciposminPreviousLine)) or
	((ciposminCurrentLine <= ciposmaxPreviousLine) and (ciposmaxCurrentLine >= ciposmaxPreviousLine)) or
	((ciposminCurrentLine >= ciposminPreviousLine) and (ciposmaxCurrentLine <= ciposmaxPreviousLine)) or
	((ciposminCurrentLine <= ciposminPreviousLine) and (ciposmaxCurrentLine >= ciposmaxPreviousLine))) or
	(svLenCurrentLine == svLenPreviousLine)):
	 	return True;

def conditionsForComparions(currentLine, previousLine):
	ciposminCurrentLine = currentLine["ciposmin"]
	ciposmaxCurrentLine = currentLine["ciposmax"]
	ciposminPreviousLine = previousLine["ciposmin"]
	ciposmaxPreviousLine = previousLine["ciposmax"]
	infoCurrentLine = currentLine["INFO"]
	infoPreviousLine = previousLine["INFO"]
	ciEndminCurrentLine = currentLine["ciendmin"]
	ciEndmaxCurrentLine = currentLine["ciendmax"]
	ciEndminPreviousLine = previousLine["ciendmin"]
	ciEndmaxPreviousLine = previousLine["ciendmax"]
	posCur = currentLine["POS"]
	posPrev = previousLine["POS"]
	endCur = infoCurrentLine["END"]
	endPrev = infoPreviousLine["END"]

	if ((infoCurrentLine["SVTYPE"] == "TRA") or (infoCurrentLine["SVTYPE"] == "BND")): #if BND or TRA check for same translocation position second chromosome
		print "tra or bnd"
		if ((infoPreviousLine["SVTYPE"] == "BND") or (infoPreviousLine["SVTYPE"] == "BND")):
			print "BND or TRA"
	 		if infoCurrentLine["SVTYPE"] == "BND":
				altCur = currentLine["ALT"]
				if "[" in altCur:
					split_string_record_alt = altCur.split("[")
				if "]" in altCur:
					split_string_record_alt = altCur.split("]")
				splitted_alt = split_string_record_alt[1].split(":")
				chrom2cur = splitted_alt[0]
				pos2cur = splitted_alt[1]
				pos2cur = int(pos2cur)
				print chrom2cur, pos2cur

			elif infoCurrentLine["SVTYPE"] == "TRA":
				chrom2cur = infoCurrentLine["chr2"]
				pos2cur = infoCurrentLine["END"]
				print chrom2cur, pos2cur

			if infoPreviousLine["SVTYPE"] == "BND":
				altPrev = previousLine["ALT"]
				if "[" in string_record_alt:
					split_string_record_alt = altPrev.split("[")
				if "]" in string_record_alt:
					split_string_record_alt = altPrev.split("]")
				splitted_alt = split_string_record_alt[1].split(":")
				chrom2prev = splitted_alt[0]
				pos2prev = splitted_alt[1]
				pos2prev = int(pos2prev)
				print chrom2prev, pos2prev

			elif infoPreviousLine["SVTYPE"] == "TRA":
				chrom2prev = infoPreviousLine["chr2"]
				pos2prev = infoPreviousLine["END"]
				print chrom2prev, pos2prev

			chrom2curMin = chrom2cur - 20
			chrom2curMax = chrom2cur + 20
			chrom2prevMin = chrom2prev - 20
			chrom2prevMax = chrom2prev + 20
			if ((chrom2cur == chrom2prev) and
			((chrom2curMax >= chrom2prevMin) and (chrom2curMin <= chrom2prevMin)) or
			((chrom2curMin <= chrom2prevMax) and (chrom2curMax >= chrom2prevMax)) or
			((chrom2curMin >= chrom2prevMin) and (chrom2curMax <= chrom2prevMax)) or
			((chrom2curMin <= chrom2prevMin) and (chrom2curMax >= chrom2prevMax)) or
			(chrom2cur == chrom2prev)):
				return True
		else:
			 print infoCurrentLine["SVTYPE"], posCur, infoPreviousLine["SVTYPE"], posPrev

	elif ((infoCurrentLine["SVTYPE"] == infoPreviousLine["SVTYPE"]) and #if types are equal
	 (currentLine["CHROM"] == previousLine["CHROM"]) and #if on the same chrom
	 (((ciposmaxCurrentLine >= ciposminPreviousLine) and (ciposminCurrentLine <= ciposminPreviousLine)) or
	 ((ciposminCurrentLine <= ciposmaxPreviousLine) and (ciposmaxCurrentLine >= ciposmaxPreviousLine)) or
	 ((ciposminCurrentLine >= ciposminPreviousLine) and (ciposmaxCurrentLine <= ciposmaxPreviousLine)) or
	 ((ciposminCurrentLine <= ciposminPreviousLine) and (ciposmaxCurrentLine >= ciposmaxPreviousLine)) or
	 (posCur == posPrev)) and #if in the same confidenceinterval of pos
	 (((ciEndmaxCurrentLine >= ciEndminPreviousLine) and (ciEndminCurrentLine <= ciEndminPreviousLine)) or
	 ((ciEndminCurrentLine <= ciEndmaxPreviousLine) and (ciEndmaxCurrentLine >= ciEndmaxPreviousLine)) or
	 ((ciEndminCurrentLine >= ciEndminPreviousLine) and (ciEndmaxCurrentLine <= ciEndmaxPreviousLine)) or
	 ((ciEndminCurrentLine <= ciEndminPreviousLine) and (ciEndmaxCurrentLine >= ciEndmaxPreviousLine)) or
	 (endCur == endPrev))): #if in the same confidencinterval of end
	 	#print infoCurrentLine["SVTYPE"], infoPreviousLine["SVTYPE"]
	 	return True

def combineVCFs(delly, manta, output):
	vcf_reader_template = vcf.Reader(filename='template.vcf')
	vcf_output_file = open(output, 'w')
	try:
		vcf_writer = vcf.Writer(vcf_output_file, vcf_reader_template)
	except IOError:
		sys.exit('Error: Cannot open vcf-file: {0}'.format(vcf_output_file))
	else:
		try:
			vcf_delly = open(delly, 'r')
			vcfD = vcf.Reader(vcf_delly)
		except IOError:
			sys.exit('Error: Cannot open vcf-file: {0}'.format(delly))

		else:
			try:
				vcf_manta = open(manta, 'r')
				vcfM = vcf.Reader(vcf_manta)
			except IOError:
				sys.exit('Error: Cannot open vcf-file: {0}'.format(manta))
			else:
				for record in vcfD:
					new_record = extractInfoDelly(record)
					vcf_writer.write_record(new_record)

				for record in vcfM:
					new_record = extractInfoManta(record)
					vcf_writer.write_record(new_record)


			vcf_delly.close()
		vcf_manta.close()
	vcf_output_file.close()

def combineVCFsOverlapFilter(delly, manta, output):
	vcf_reader_template = vcf.Reader(filename='template.vcf')
	vcf_output_file = open(output, 'w')
	try:
		vcf_writer = vcf.Writer(vcf_output_file, vcf_reader_template)
	except IOError:
		sys.exit('Error: Cannot open vcf-file: {0}'.format(vcf_output_file))
	else:
		try:
			vcf_delly = open(delly, 'r')
			vcfD = vcf.Reader(vcf_delly)
		except IOError:
			sys.exit('Error: Cannot open vcf-file: {0}'.format(delly))

		else:
			try:
				vcf_manta = open(manta, 'r')
				vcfM = vcf.Reader(vcf_manta)
			except IOError:
				sys.exit('Error: Cannot open vcf-file: {0}'.format(manta))
			else:
				list_all_records_MantaDelly = []
				for record in vcfD:
					new_recordDelly = extractInfoDelly(record)
					list_Delly = getInfoInList(new_recordDelly)
					list_all_records_MantaDelly.append(list_Delly)

				for record in vcfM:
					new_recordManta = extractInfoManta(record)
					list_Manta = getInfoInList(new_recordManta)
					list_all_records_MantaDelly.append(list_Manta)

				#first we will sort on the chrom, then on the startposition of the SV
				list_all_records_MantaDelly.sort(key=operator.itemgetter("CHROM","POS"))
				list_for_comparisonSVs = []
				all_svs_to_write = []

				previousLine = list_all_records_MantaDelly[0]
				for currentLine in list_all_records_MantaDelly:
					if (conditionsInsForComparison(currentLine, previousLine)):
						list_for_comparisonSVs.append(currentLine)

					elif (conditionsForComparions(currentLine, previousLine)):
						list_for_comparisonSVs.append(currentLine)

					else:
						list_sv_to_print = compareFilterSVs(list_for_comparisonSVs)
						all_svs_to_write.extend(list_sv_to_print)
						list_for_comparisonSVs = startNewListForComparisonSVs(list_for_comparisonSVs, currentLine)

					previousLine = currentLine

				if len(list_for_comparisonSVs) == 1: #when the last line in the file was a single SV
					list_sv_to_print = compareFilterSVs(list_for_comparisonSVs)
					all_svs_to_write.extend(list_sv_to_print)

				filterBndAndWriteSVs(all_svs_to_write, vcf_writer)

			vcf_delly.close()
		vcf_manta.close()
	vcf_output_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Merging manta and delly file of same BAM')

	parser.add_argument('--filterOverlap', help = "outputVCF contains SVs only once. When SV is called twice, only manta output is given. ", action = "store_true" )

	required_named = parser.add_argument_group('Required arguments')
	required_named.add_argument('-d', '--dellyVCF', help = "input a vcf file from delly to process", required=True)
	required_named.add_argument('-m', '--mantaVCF', help = "input a vcf file from manta to process", required=True)
	required_named.add_argument('-o', '--outputVCF', help = "give the filename of a vcf file for the merged output", required=True)

	args = parser.parse_args()
	if args.filterOverlap:
		print "overlap function used, no double values printed"
		combineVCFsOverlapFilter(args.dellyVCF, args.mantaVCF, args.outputVCF)
	else:
		combineVCFs(args.dellyVCF, args.mantaVCF, args.outputVCF)

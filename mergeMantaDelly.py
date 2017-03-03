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
	elif record.INFO["SVTYPE"] == "BND": #because BNDs do not have a value for END and contain mostly TRA with a length of 0, END is set equal to POS.
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

	if "CIPOS" in new_record.INFO:
		cipos = new_record.INFO["CIPOS"]
		ciposmin = pos + cipos[0] - 20
		ciposmax = pos + cipos[1] + 20
	else:
		cipos = [-50, 50]
		ciposmin = pos + cipos[0]
		ciposmax = pos + cipos[1]

	if "CIEND" in new_record.INFO:
		ciend= new_record.INFO["CIEND"]
		ciendmin = (end + ciend[0] - 20)
		ciendmax = (end + ciend[1] + 20)
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

	dict_record = {"ciposint": ciposint, "tool": tool, "CHROM" : chrom, "POS": pos, "ciposmin": ciposmin, "ciposmax": ciposmax, "ciendmin": ciendmin, "ciendmax": ciendmax, "ID" : ID, "INFO": info, "record": new_record}
	return dict_record

def startNewListForComparisonSVs(list_for_comparisonSVs, currentLine):
	comparingSVs_list = list_for_comparisonSVs
	if comparingSVs_list:
		comparingSVs_list = []
	comparingSVs_list.append(currentLine)

	return comparingSVs_list

def compareFilterAndWriteSVs(list_for_comparisonSVs, vcf_writer): #is called when list_for_comparisonSVs is about to be refreshed
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
			if similar_SVs[1]["tool"] == "DELLY":
				record_sv_to_print.INFO["CSA"] = 2
				record_sv_to_print.INFO["INFODELLY"] = similar_SVs[1]["INFO"]
			else:
				record_sv_to_print.INFO["CSA"] = 1

			sv_to_print.append(record_sv_to_print)
		else:
			record_sv_to_print = similar_SVs[0]["record"]
			record_sv_to_print.INFO["CSA"] = 1
			sv_to_print.append(record_sv_to_print)

	for record in sv_to_print:
		vcf_writer.write_record(record)

def conditionsInsForComparison(currentLine, previousLine):
	#for insertions there is checked on length instead of END of SV
	svLenCurrentLine = currentLine["INFO"]["SVLEN"]
	svLenPreviousLine = previousLine["INFO"]["SVLEN"]
	cilength = [-20,20]
	svtypeCur = currentLine["INFO"]["SVTYPE"]
	svtypePrev = previousLine["INFO"]["SVTYPE"]

	if (((svtypeCur == "INS") and (svtypePrev == "INS")) and #for insertions
	(((svLenCurrentLine + cilength[0]) <= (svLenPreviousLine + cilength[1])) or #if lengths are in the same confidenceinterval
	((svLenCurrentLine + cilength[1]) >= (svLenPreviousLine + cilength[0])))):
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

	if (((infoCurrentLine["SVTYPE"] == infoPreviousLine["SVTYPE"]) or (infoCurrentLine["SVTYPE"] == "TRA" and infoPreviousLine["SVTYPE"] == "BND") or
	 (infoCurrentLine["SVTYPE"] == "BND" and infoPreviousLine["SVTYPE"] == "TRA")) and ((currentLine["CHROM"] == previousLine["CHROM"]) and #if on the same chrom
	 (ciposmaxCurrentLine >= ciposminPreviousLine or ciposminCurrentLine <= ciposmaxPreviousLine) and  #if positions are in the same confidence interval
	 (ciEndmaxCurrentLine >= ciEndminPreviousLine or ciEndminCurrentLine <= ciEndmaxPreviousLine))): #if ends are in the same confidence interval
	 	return True;


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

				previousLine = list_all_records_MantaDelly[0]
				for currentLine in list_all_records_MantaDelly:
					if (conditionsInsForComparison(currentLine, previousLine)):
						list_for_comparisonSVs.append(currentLine)

					elif (conditionsForComparions(currentLine, previousLine)):
						list_for_comparisonSVs.append(currentLine)

					else:
						compareFilterAndWriteSVs(list_for_comparisonSVs, vcf_writer)
						list_for_comparisonSVs = startNewListForComparisonSVs(list_for_comparisonSVs, currentLine)

					previousLine = currentLine

				if len(list_for_comparisonSVs) == 1: #when the last line in the file was a single SV
					compareFilterAndWriteSVs(list_for_comparisonSVs, vcf_writer)

			vcf_delly.close()
		vcf_manta.close()
	vcf_output_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Merging manta and delly file of same BAM')

	parser.add_argument('--filterOverlap', help = "outputVCF contains SVs only once. When SV is called twice, only manta output is given. ", action = "store_true" )

	required_named = parser.add_argument_group('Required arguments')
	required_named.add_argument('-d', '--dellyVCF', help = "input a vcf file from delly to process", required=True)
	required_named.add_argument('-ma', '--mantaVCF', help = "input a vcf file from manta to process", required=True)
	required_named.add_argument('-o', '--outputVCF', help = "give the filename of a vcf file for the merged output", required=True)

	args = parser.parse_args()
	if args.filterOverlap:
		print "overlap function used, no double values printed"
		combineVCFsOverlapFilter(args.dellyVCF, args.mantaVCF, args.outputVCF)
	else:
		combineVCFs(args.dellyVCF, args.mantaVCF, args.outputVCF)

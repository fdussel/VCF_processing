#template.vcf is needed to make a merged file
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
	new_record.INFO["END"] = end

	length = 0

	if "SVLEN" in record.INFO:
		length = record.INFO["SVLEN"]
	elif record.INFO["SVTYPE"] == "BND":
		length = 0
	elif record.INFO["SVTYPE"] == "INS":
		length = 0
		print record.INFO["SVTYPE"], record.CHROM, record.POS, "-length of insertion unknown, set to zero"
	new_record.INFO["SVLEN"] = length

	return new_record


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
				chromManta = 0
				posManta = 0
				list_Manta_per_record = []

				for record in vcfM:
					recordManta = record
					new_recordManta = extractInfoManta(record)
					vcf_writer.write_record(new_recordManta)

					chromManta = record.CHROM
					posManta = record.POS
					#confidence interval position
					if "CIPOS" in record.INFO:
						ciposManta = record.INFO["CIPOS"]
					else:
						ciposManta = [-50, 50]
					ciposMantaMin = posManta + ciposManta[0]
					ciposMantaMax = posManta + ciposManta[1]
					#dict_manta_per_record = {"chromManta": chromManta, "ciposMantaMin": ciposMantaMin, "ciposMantaMax": ciposMantaMax}
					list_record_info = [chromManta, ciposMantaMin, ciposMantaMax]
					list_Manta_per_record.append(list_record_info)

				for record in vcfD:
					recordDelly = record
					chromDelly = record.CHROM
					posDelly = record.POS
					ciposDelly = record.INFO["CIPOS"]
					ciposDellyMin = posDelly + ciposDelly[0]
					ciposDellyMax = posDelly + ciposDelly[1]
					new_recordDelly = extractInfoDelly(recordDelly)

					#ciposMantaMin = dict_manta_per_record["ciposMantaMin"]
					#ciposMantaMax = dict_manta_per_record["ciposMantaMax"]
					#chromManta = dict_manta_per_record["chromManta"]
					recordwritten = False
					#print list_Manta_per_record
					for list_record_info in list_Manta_per_record:
						#print list_record_info
						chromManta = list_record_info[0]
						ciposMantaMin = list_record_info[1]
						ciposMantaMax = list_record_info[2]
						#print chromDelly
						#if chromDelly == chromManta:
						#	print ""

						 	#and recordwritten == False:
							#vcf_writer.write_record(new_recordDelly)
							#recordwritten = True
						#else:
						#	if (ciposDellyMax <= ciposMantaMin or ciposDellyMin >= ciposMantaMax) and recordwritten == False:
						#		vcf_writer.write_record(new_recordDelly)
						#		recordwritten = True

			vcf_delly.close()
		vcf_manta.close()
	vcf_output_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Merging manta and delly file of same BAM')

	parser.add_argument('--overlapFilter', help = "outputVCF contains SVs only once. When SV is called twice, only manta output is given. ", action = "store_true" )

	required_named = parser.add_argument_group('Required arguments')
	required_named.add_argument('-d', '--dellyVCF', help = "input a vcf file from delly to process", required=True)
	required_named.add_argument('-ma', '--mantaVCF', help = "input a vcf file from manta to process", required=True)
	required_named.add_argument('-o', '--outputVCF', help = "give the filename of a vcf file for the merged output", required=True)

	args = parser.parse_args()
	if args.overlapFilter:
		print "overlap function used, no double values printed"
		combineVCFsOverlapFilter(args.dellyVCF, args.mantaVCF, args.outputVCF)
	else:
		combineVCFs(args.dellyVCF, args.mantaVCF, args.outputVCF)

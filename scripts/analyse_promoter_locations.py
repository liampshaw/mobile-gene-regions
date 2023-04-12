import argparse
import pandas as pd
import re

def get_options():
    parser = argparse.ArgumentParser(description="Analyse promoter locations relative to a region of a genome",
                                     prog="analyse_promoter_locations")
    parser.add_argument("--genome", help="input genome", required=True)
    parser.add_argument("--flanking_size", help="size of flanking region upstream/downstream",
    	required=False, default=5000)
    parser.add_argument("--central_region", help="central gene region to use (in coordinates of whole genome). format: x-y or complement(x-y)", required=False, default="none")
    parser.add_argument("--fasta_dir", help="directory where genome is stored in genome.fa", required=True)
    parser.add_argument("--promoter_dir", help="directory where promoter information is stored in genome-promoter_predictions.tsv", required=True)
    parser.add_argument("--output", help="output file (csv)", required=False, default="")
    #parser.add_argument("--reverse_complement", "whether to reverse complement the results", action="store_true", default=False, required=False)
    return parser.parse_args()


def extract_region_from_df(df, start_limit, end_limit):
	"""extracts all rows of df with centre of promoter within limits """
	cut_one = df[df["location"]>start_limit]
	cut_two = cut_one[cut_one["location"]<end_limit]
	return(cut_two)

# NZ_AP021952.1 864bp 27208-28072

def main():
	args = get_options()
	# by default, assume not complement and positive strand unless we hear otherwise
	complement = False
	strand = "+"

	genome = args.genome
	promoter_dir = args.promoter_dir

	# The information from the genomes needed is typically stored in output/analysis/sequence/{focal_gene}_seqs_extracted_from_contigs.fa
	promoter_locations_df = pd.read_csv(promoter_dir+"/"+genome+"-promoter_predictions.tsv",
		sep="\t", header=0)

	central_region = args.central_region
	flanking_size = int(args.flanking_size)

	if "complement" in central_region:
		complement = True
		strand = "-"

	if complement==True:
		central_region_limits = [int(x) for x in re.sub("[()]", "", re.sub("complement", "", central_region)).split("-")[::-1]]
	elif complement==False:
		central_region_limits = [int(x) for x in central_region.split("-")]
	# bear in mind: need to handle reverse complement

	# For handling the strand, we select only those with the right strand
	# Promoters going in wrong direction we are not interested in
	promoter_locations_df = promoter_locations_df[promoter_locations_df["strand"]==strand]
	# add average location of promoter (centre)
	# Seems better that I use the centre of the promoter rather than the start/end
	# Promoter is 40nt stretch. In the df it is always stored start<end independent of strand
	# Otherwise introduce a systematic bias upstream/downstream (a small one by 20nt)
	promoter_locations_df = promoter_locations_df.assign(location= (promoter_locations_df["start"]+20))#promoter_locations_df["end"])/2)

	# Calculate the requested limits
	# and then extract the regions that satisfy
	if complement==False:
		upstream_region = [central_region_limits[0]-flanking_size, central_region_limits[0]]
		downstream_region = [central_region_limits[1], central_region_limits[1]+flanking_size]
	elif complement==True:
		upstream_region = [central_region_limits[0], central_region_limits[0]+flanking_size]
		downstream_region = [central_region_limits[1]-flanking_size, central_region_limits[1]]

	downstream_promoters_df = extract_region_from_df(promoter_locations_df, downstream_region[0], downstream_region[1])
	upstream_promoters_df = extract_region_from_df(promoter_locations_df, upstream_region[0], upstream_region[1])

	# Measure the distance between the 'end' of the promoter sequence and the gene start/end
	if complement==False:
		downstream_promoters_df = downstream_promoters_df.assign(distance=downstream_promoters_df["location"]-central_region_limits[1])
		upstream_promoters_df = upstream_promoters_df.assign(distance=-(central_region_limits[0]-upstream_promoters_df["location"]))
	elif complement==True:
		# the strand is negative, so the 'end' of the promoter is really its start and vice versa
		downstream_promoters_df = downstream_promoters_df.assign(distance=central_region_limits[1]-downstream_promoters_df["location"])
		upstream_promoters_df = upstream_promoters_df.assign(distance=-(upstream_promoters_df["location"]-central_region_limits[0]))

	output_df = pd.concat([downstream_promoters_df,upstream_promoters_df])
	if args.output=="":
		print(output_df)
	else:
		output_df.to_csv(args.output, index=False)

# Q: is there a risk that I miss promoters that overlap with the start of 
# the gene itself by how I'm selecting these?

# Now what do we do with these data frames?

# Can 

if __name__=="__main__":
	main()
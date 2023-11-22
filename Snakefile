#configurations for running the code
#Anything you would normally "hardcode" like paths, variables for the dada2 pipeline, etc go here
configfile: "config.yaml"

#Gets sample names from the forward reads
WC = glob_wildcards(config['path']+"{names}_1.fastq.gz") # CHANGE THIS
SAMPLES = WC.names
R1 = expand(config['path']+'{names}_1.fastq.gz', names=SAMPLES) #CHANGE THIS
R2 = expand(config['path']+'{names}_2.fastq.gz', names=SAMPLES) #CHANGE THIS

#local rules marks a rule as local and does not need to be submitted as a job to the cluster
localrules: all, plotQP, track

#this rule specifies all things you want generated
rule all:
	input:
		"output/errorRates_R1.rds",
		"output/seqtab_nochimeras.rds",
		"output/ASVs.txt",
		"output/taxonomy.txt",
		"output/ASVseqs.txt",
		"output/tree.nwk",
		"output/track.tsv"

#clears all outputs (except for plotted quality profiles)
rule clean:
    shell:
        '''
        rm output/*
	rm logs/*
        '''

#plots quality profiles
rule plotQP:
	input: R1, R2
	output:
		R1 = expand('output/qualityProfiles/R1/{sample}_R1_qual.jpg',sample=SAMPLES),
		R2 = expand('output/qualityProfiles/R2/{sample}_R2_qual.jpg',sample=SAMPLES)
	script:
		'scripts/plotQP.R'

#quality filters R1 and R2 (forward and reverse reads)
rule filter:
	input:
		R1=R1,
		R2=R2
	output:
		R1 = expand(config['path']+"filtered/{sample}_R1.fastq.gz", sample=SAMPLES),
		R2 = expand(config['path']+"filtered/{sample}_R2.fastq.gz", sample=SAMPLES),
		track = temp("output/track_filtered.txt")
	params:
		samples = SAMPLES
	log:
		"logs/filter.txt" #I always have the logs go to one place so I can easily see what went wrong
	script:
		"scripts/filter.R"

#error modeling and plotting the errors
rule learnError:
	input:
		R1 = rules.filter.output.R1,
		R2 = rules.filter.output.R2 #note you can declare the output of another rule as a dependency
	output:
		errR1 = "output/errorRates_R1.rds",
		errR2 = "output/errorRates_R2.rds",
		plotErrR1 = "output/figures/errorRates_R1.pdf",
		plotErrR2 = "output/figures/errorRates_R2.pdf"
	log:
		"logs/learnError.txt"
	script:
		"scripts/learnError.R"

#dereplicates and merges forward and reverse reads, finding all unique sequences
rule dereplicate:
	input:
		R1 = rules.filter.output.R1,
		R2 = rules.filter.output.R2,
		errR1 = rules.learnError.output.errR1,
		errR2 = rules.learnError.output.errR2
	output:
		seqtab = "output/seqtab_withchimeras.rds",
		track = temp("output/track_derep.txt")
	params:
		samples = WC.names
	log:
		"logs/dereplicate.txt"
	script:
		"scripts/dereplicate.R"

#this is where the chimeras get removed
rule removeChimeras:
	input:
		seqtab = rules.dereplicate.output.seqtab,
	output:
		seqtab = "output/seqtab_nochimeras.rds",
		track = temp("output/track_removeChimeras.txt")
	params:
		samples = WC.names
	log:
		"logs/removeChimeras.txt"
	script:
		"scripts/removeChimeras.R"

#this is how the reads are tracked
rule track:
	input:
		filt = rules.filter.output.track,
		derep = rules.dereplicate.output.track,
		nochim = rules.removeChimeras.output.track
	output:
		track = "output/track.tsv"
	log:
		"logs/track.txt"
	script:
		"scripts/track.R"

#assignes taxonomy using the silva database
rule taxonomy:
	input:
		seqtab = rules.removeChimeras.output.seqtab,
	output:
		otus = "output/ASVs.txt",
		taxonomy = "output/taxonomy.txt",
		ASVseqs = "output/ASVseqs.txt"
	params:
		samples = SAMPLES
	log:
		"logs/taxonomy.txt"
	script:
		"scripts/taxonomy.R"

#Generate de novo phylogenetic trees
rule alignment:
	input:
		seqtab = rules.removeChimeras.output.seqtab
	output:
		alignment = "output/alignment.fasta",
	log:
		"logs/alignment.txt"
	script:
		"scripts/alignment.R"

rule tree:
	input:
		alignment = rules.alignment.output.alignment
	output:
		tree = "output/tree.nwk"
	log:
		"logs/tree.txt"
	shell:
		"""
		set +u
		module load anaconda/5.1
		source activate snakemakedada
		set -u
		FastTreeMP -gamma -nt -log {log} < {input} > {output}
		"""

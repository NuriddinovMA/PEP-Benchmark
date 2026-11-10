from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers

import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="the script to predict mutation effect on the gene expression")
parser.add_argument('--regulome','-r', default="peps/EP.GRCh38-c.tsv", help="he path to files with p-e pairs")
parser.add_argument('--ref_genome','-g', default="genomes/hg38.fa", help="the path to a genome file")
parser.add_argument('--out','-o', default="alpha_score", help="the file name for the output")
parser.add_argument('--batch','-b', type=int, default=100, help="the number of sequence in batch")
parser.add_argument('--threads','-t', type=int, default=5, help="the threads count")
parser.add_argument('--step','-s', type=int, default=0, help="from which PE-pair must start the prediction. PEP = batch*step")
parser.add_argument('--end','-e', type=int, default=False, help="on which PE-pair must end the prediction. PEP = batch*end")
parser.add_argument('--window','-w', choices=[1048576, 524288, 131072,16384], type=int, help="the length of sequence to prediction")
parser.add_argument('--use_dropped','-u', action='store_true', help="if true, use the peps even the enhancer don't fall in window")
args=parser.parse_args()

API_KEY = 'AIzaSyBiTXcVfXvNx5N1Fz0TF3_8XyiWkAKI2AQ'
model = dna_client.create(API_KEY)
window = args.window // 2
batch = args.batch
FA = SeqIO.to_dict(SeqIO.parse(args.ref_genome, 'fasta'))

eps_hash = {}
with open(args.regulome,'r') as f: lines = f.readlines()[1:]
for line in lines:
	parse = line.split()
	chrm_ep, start_ep, end_ep, mut_start, gene_id, chr_tss, start_tss, end_tss = parse[:8]
	start_ep, end_ep, mut_start, start_tss, end_tss = int(start_ep), int(end_ep), int(mut_start),int(start_tss), int(end_tss)
	if parse[12] == "TRUE": regulate = 1
	else: regulate = 0
	if args.use_dropped or ((start_tss-window < mut_start) and (start_tss+window > mut_start+2000)):
		eps_hash[gene_id,start_tss,start_ep] = ( chr_tss, start_tss, chrm_ep, mut_start, regulate)

gene_ids = sorted(eps_hash)
if args.end: pass
else: args.end = len(gene_ids)//batch + 2
print(len(gene_ids))

for step in range(args.step, args.end):
	print(step,file=sys.stderr)
	all_df_dropped = []
	all_df_scored = []
	droped = {}
	
	variants,intervals,scored_intervals = [],[],[]
	for gene_id,start_tss,start_ep in gene_ids[(step*batch):((step+1)*batch)]:
		chr_tss, start_tss, chrm_ep, mut_start, regulate = eps_hash[gene_id,start_tss,start_ep ]
		if (mut_start > start_tss-window) and (start_tss+window > mut_start+2000):
			scored_interval = f"{chr_tss}:{start_tss-window}-{start_tss+window}:."
			scored_intervals.append( (gene_id, scored_interval) )
			interval = genome.Interval(chromosome=chr_tss, start=start_tss-window, end=start_tss+window)
			print(f"{chrm_ep}:{mut_start}-{mut_start+2000} {chr_tss}:{start_tss-window}-{start_tss+window}",file=sys.stderr)
			reference = str(FA[chrm_ep].seq[(mut_start):(mut_start+2000)]).upper()
			shuffled = np.array([list(reference)]*10)
			np.apply_along_axis(np.random.shuffle,1,shuffled)
			shuffled = np.apply_along_axis(''.join,1,shuffled).tolist()
			scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
			for j in range(10):
				variant = genome.Variant(
					chromosome=chrm_ep,
					position=mut_start,
					reference_bases=reference,
					alternate_bases=shuffled[j],
					name=f"{chrm_ep}:{start_ep}:{end_ep}",
				)
				variants.append(variant)
				intervals.append(interval)
		elif args.use_dropped:
			droped['gene_id'].append(gene_id)
			droped['variant_id'].append(f"{chrm_ep}:{mut_start}")
			droped['scored_interval'].append(scored_interval)
			droped['raw_score'].append(0)
			droped['quantile_score'].append(0)
		else: pass
	
	variant_scores = model.score_variants(
		intervals=intervals,
		variants=variants,
		variant_scorers=[scorer],
		progress_bar = False,
		max_workers=args.threads,
	) 
	try:
		df_scored = variant_scorers.tidy_scores(variant_scores)
		all_df_scored.append(df_scored)
	except TypeError: pass
	
	if args.use_dropped:
		all_df_dropped.append(pd.DataFrame.from_dict(droped))
		try: df_dropped = pd.concat(all_df_dropped)
		except ValueError: pass
		
	try:
		df_scored = pd.concat(all_df_scored)
		df_scored['scored_interval'] = df_scored['scored_interval'].astype('string')
		a = df_scored[['gene_id','scored_interval']].apply(tuple,1).isin(scored_intervals)
		df_scored = df_scored[a]
		df_short = df_scored[['gene_id','variant_id','scored_interval','raw_score','quantile_score']][(df_scored['Assay title'] == "total RNA-seq") & (df_scored['ontology_curie'] == "EFO:0002067")]
		df_short['variant_id'] = df_short['variant_id'].astype('string').str.split(':')
		df_short['variant_id'] = df_short['variant_id'].apply(lambda x: f"{x[0]}:{int(x[1])}")
	except ValueError: pass
	
	e1,e2 = True,True
	try: df_short
	except NameError: e1 = False
	try: df_dropped
	except NameError: e2 = False
	
	if e1 and e2: df = pd.concat([df_short,df_dropped])
	elif e1: df = df_short
	elif e2: df = df_dropped
	else: break
	
	df.to_csv(f'scores/{args.out}.{step}.tsv',sep='\t',index=False)
	del df_scored, df_short, df

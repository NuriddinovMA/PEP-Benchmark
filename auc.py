import glob
import sys
import argparse
import numpy as np
import sklearn.metrics as skm

parser = argparse.ArgumentParser(description="the script to predict mutation effect on the gene expression")
parser.add_argument('--regulome','-r', default="peps/EP.GRCh38-c.tsv", help="the path to files with p-e pairs")
parser.add_argument('--path','-p', default="scores", help="the file name for the output")
parser.add_argument('--name','-n', default="alpha2_score", help="the base file name for generated predictionss")
parser.add_argument('--min_dists','-m', nargs='*', type=int,default=[0], help="the base file name for generated predictionss")
parser.add_argument('--max_dists','-M', nargs='*', type=int,default=[2500000], help="the base file name for generated predictionss")
parser.add_argument('--out','-o', default="alpha2_auc", help="the file name for the output")
args=parser.parse_args()

train = args.regulome
predictions = args.name
min_dists = args.min_dists
max_dists = args.max_dists
ln_dists = len(min_dists)
assert len(min_dists) == len(max_dists)
dists_list = [(min_dists[i],max_dists[i]) for i in range(ln_dists)]

dists_list = [(min_dists[i],max_dists[i]) for i in range(ln_dists)]
Train,Dists = {},{}
with open(train,'r') as f: lines = f.readlines()
for line in lines[1:]:
	parse = line.split()
	chrm_ep, start_ep, end_ep, mut_start, gene_id, chr_tss, start_tss, end_tss = parse[:8]
	start_ep, end_ep, mut_start, start_tss, end_tss = int(start_ep), int(end_ep), int(mut_start),int(start_tss), int(end_tss)
	dist = int(parse[9])
	if parse[11] == "TRUE": regulate = 1
	else: regulate = 0
	try: 
		Train[gene_id][chrm_ep, mut_start] = regulate
		Dists[gene_id][chrm_ep, mut_start] = dist
	except KeyError: 
		Train[gene_id] = {(chrm_ep, mut_start): regulate }
		Dists[gene_id] = {(chrm_ep, mut_start): dist }


Pred = {}
all_raw_score, all_perc_score = [],[]
files = glob.glob(f'{args.path}/{args.name}.*.tsv')
for file in files:
	with open(f"{file}",'r') as f: lines = f.readlines()
	for line in lines[1:]:
		gene_id, chrm_center, chrm_locus, raw_score, perc_score = line.split()
		raw_score, perc_score = float(raw_score), float(perc_score)
		chrm_ep,center = chrm_center.split(':')
		center = int(center)
		try: 
			Pred[gene_id]
			try: 
				Pred[gene_id][chrm_ep,center][0].append(raw_score)
				Pred[gene_id][chrm_ep,center][1].append(perc_score)
			except KeyError: Pred[gene_id][chrm_ep,center] = [[raw_score],[perc_score]]
		except KeyError: Pred[gene_id] = {(chrm_ep,center): [[raw_score],[perc_score]]}
		
for gene_id in Pred:
	for chrm_ep,center in Pred[gene_id]:
		Pred[gene_id][chrm_ep,center][0] = np.mean( np.abs(Pred[gene_id][chrm_ep,center][0]) )
		Pred[gene_id][chrm_ep,center][1] = np.mean( np.abs(Pred[gene_id][chrm_ep,center][1]) )
		all_raw_score.append(Pred[gene_id][chrm_ep,center][0])
		all_perc_score.append(Pred[gene_id][chrm_ep,center][1])

TP_all = []
FN_all = []
FP_all = []

out=open(f"{args.path}/{args.out}.tsv",'w')
out.write(f"metric\tdistances\traw_variant_scoring\tperc_variant_scoring\n")
reg = { (min_dists[i],max_dists[i]):[] for i in range(ln_dists) }
score_raw = { (min_dists[i],max_dists[i]):[] for i in range(ln_dists) }
score_perc = { (min_dists[i],max_dists[i]):[] for i in range(ln_dists) }

for gene_id in Train:
	for chrm_ep,center in Train[gene_id]:
		regulate = Train[gene_id][chrm_ep,center]
		dist = Dists[gene_id][chrm_ep,center]
		for i,(min_dist,max_dist) in enumerate(dists_list):
			if min_dist <= dist <= max_dist:
				try: 
					Pred[gene_id]
					try: 
						raw_effect,perc_effect = Pred[gene_id][chrm_ep,center]
						score_raw[min_dist,max_dist].append(raw_effect)
						score_perc[min_dist,max_dist].append(perc_effect)
						reg[min_dist,max_dist].append(regulate)
					except KeyError: pass#print(gene_id,Pred[gene_id],chrm_ep,center,file=sys.stderr)
				except KeyError: pass#print(gene_id,chrm_ep,center)
"precn = TP / (TP + FP)"
"recll = TP / (TP + FN)"

for i,(min_dist,max_dist) in enumerate(dists_list):
	precision_raw, recall_raw, thresholds_raw = skm.precision_recall_curve(reg[min_dist,max_dist],score_raw[min_dist,max_dist])
	precision_perc, recall_perc, thresholds_perc = skm.precision_recall_curve(reg[min_dist,max_dist],score_perc[min_dist,max_dist])
	auc_raw, auc_perc = skm.auc(recall_raw,precision_raw),skm.auc(recall_perc,precision_perc)
	out.write(f"auc\t{min_dist}-{max_dist}\t{auc_raw}\t{auc_perc}\n")
out.close()
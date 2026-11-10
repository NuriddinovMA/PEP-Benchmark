import glob
import sys
import argparse
import sklearn.metrics as skm


parser = argparse.ArgumentParser(description="the script to predict mutation effect on the gene expression")
parser.add_argument('--regulome','-r', default="peps/EP.GRCh38-c.tsv", help="the path to files with p-e pairs")
parser.add_argument('--name','-n', default="alpha_score", help="the base file name for generated predictions")
parser.add_argument('--out','-o', default="alpha_auc", help="the file name for the output")
rgs=parser.parse_args()

train = args.regulome
prediction = args.name

Train = {}
with open(train,'r') as f: lines = f.readlines()
for line in lines[1:]:
	parse = line.split()
	chrm_ep, start_ep, end_ep, mut_start, gene_id, chr_tss, start_tss, end_tss = parse[:8]
	start_ep, end_ep, mut_start, start_tss, end_tss = int(start_ep), int(end_ep), int(mut_start),int(start_tss), int(end_tss)
	dist = int(parse[9])
	if parse[11] == "TRUE": regulate = 1
	else: regulate = 0
	if dist < 6000:
		try: Train[gene_id][chrm_ep, mut_start] = regulate
		except KeyError: Train[gene_id] = {(chrm_ep, mut_start): regulate}

Pred = {}
all_raw_score, all_perc_score = [],[]
files = glob.glob(f'scores/{prediction}.*.tsv')
for file in files:
	with open(f"scores/{file},'r') as f: lines = f.readlines()
	for line in lines[:]:
		gene_id, chrm_center, chrm_locus, raw_score, perc_score = line.split()
		# print(line, raw_score, perc_score)
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
raw_step = np.max(all_raw_score)/100.
perc_step = np.max(all_perc_score)/100.
reg,prc_raw,prc_perc = [],[],[]
posit_gene,posit_enh = [],[]
for gene_id in Train:
	for chrm_ep,center in Train[gene_id]:
		regulate = Train[gene_id][chrm_ep,center]
		try: 
			Pred[gene_id]
			try: 
				raw_effect, perc_effect = Pred[gene_id][chrm_ep,center]
				prc_raw.append(raw_effect)
				prc_perc.append(perc_effect)
				reg.append(regulate)
			except KeyError: print(gene_id,Pred[gene_id],chrm_ep,center,file=sys.stderr)
		except KeyError: print(gene_id,chrm_ep,center)

"precn = TP / (TP + FP)"
"recll = TP / (TP + FN)"

precision_raw, recall_raw, thresholds_raw = skm.precision_recall_curve(reg,prc_raw)
precision_perc, recall_perc, thresholds_perc = skm.precision_recall_curve(reg,prc_perc)
auc_raw, auc_perc = skm.auc(recall_raw,precision_raw),skm.auc(recall_perc,precision_perc)
with open(f"scores/{args.out}.tsv") as f:
	f.write(f"metric\traw_variant_scoring\tpercentile_variant_scoring}\n")
	f.write(f"presicion\t{precision_raw}\t{precision_perc}\n")
	f.write(f"recall\t{recall_raw}\t{recall_perc}\n")
	f.write(f"auc\t{auc_raw}\t{auc_perc}\n")

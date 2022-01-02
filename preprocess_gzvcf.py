'''
Author: Sumaiya Nazeen (sumaiya_nazeen@hms.harvard.edu)
This program prepares the input mutations file expected by the RVTT program from a filtered gzipped vcf file.
Input:
- filtered gzipped vcf file: The vcf file must be annotated with SnpEff, SnpSIFT and dbNSFP v4.0a database.
- tab-separated phenotype file in .fam format
- minimum minor allele frequency cutoff (suggested value = 0)
- maximum minor allele frequency cutoff (suggested value = 0.05)
- prefix of the output file
Output: 
- outprefix_mutations.tsv file that contains the necessary columns for running RVTT
'''


import sys
import numpy as np
import pandas as pd
from collections import Counter
import gzip

damaging = ['frameshift_variant','start_lost','stop_gained','stop_lost','splice_donor_variant','splice_acceptor_variant','splice_region_variant','structural_interaction_variant','initiator_codon_variant', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion']
excl = ['3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant', 'bidirectional_gene_fusion', 'downstream_gene_variant', 'gene_fusion',  'intergenic_region', 'intron_variant', 'non_coding_transcript_exon_variant', 'sequence_feature', 'TF_binding_site_variant', 'upstream_gene_variant']
incl = ['missense_variant','synonymous_variant'] + damaging
crit1 = ['mis_vs_syn','mis_vs_all','damaging_vs_all','damaging_vs_syn_nut']
crit2 = ['rare','all']

def print_matrix(inmap, indivs, genes, outfile):
	genes = list(genes)
	in_genes = list(inmap.keys())
	mat = np.zeros((len(indivs),len(genes)))
	for k in in_genes:
		c = Counter(inmap[k])
		c_indivs = list(c.keys())
		ind = [indivs.index(v) for v in c_indivs]
		for i in range(len(c_indivs)):
			mat[ind[i],genes.index(k)] = c[c_indivs[i]]
	df = pd.DataFrame(mat, index=indivs, columns=genes)
	df.to_csv(outfile, sep=',')
		
			
def israre(gnomgAF, gnomgAF_nfe, gnomeAF, gnomeAF_nfe, af, cutoff):
	print(gnomgAF, gnomgAF_nfe, gnomeAF, gnomeAF_nfe, af, cutoff)
	if gnomgAF == ".":
		gnomgAF = 0
	if gnomgAF_nfe == ".":
		gnomgAF_nfe = 0
	if gnomeAF == ".":
		gnomeAF = 0
	if gnomeAF_nfe == ".":
		gnomeAF_nfe = 0
	if gnomgAF == 0 and gnomgAF_nfe==0 and gnomeAF==0 and gnomeAF_nfe == 0:
		if af <= cutoff:
			return True
		else:
			return False
	else:
		if gnomgAF <= cutoff or gnomgAF_nfe <= cutoff or gnomeAF <= cutoff or gnomeAF_nfe <= cutoff:
			return True
		else:
			return False



def isInRange(gnomgAF, gnomgAF_nfe, gnomeAF, gnomeAF_nfe, af, mincutoff, maxcutoff):
        print(gnomgAF, gnomgAF_nfe, gnomeAF, gnomeAF_nfe, af, mincutoff, maxcutoff)
        if gnomgAF == ".":
                gnomgAF = 0
        if gnomgAF_nfe == ".":
                gnomgAF_nfe = 0
        if gnomeAF == ".":
                gnomeAF = 0
        if gnomeAF_nfe == ".":
                gnomeAF_nfe = 0
        if gnomgAF == 0 and gnomgAF_nfe==0 and gnomeAF==0 and gnomeAF_nfe == 0:
                if af >= mincutoff and af <= maxcutoff:
                        return True
                else:
                        return False
        else:
                if (gnomgAF>= mincutoff and gnomgAF <= maxcutoff) or (gnomgAF_nfe >= mincutoff and gnomgAF_nfe <= maxcutoff) or (gnomeAF >= mincutoff and gnomeAF <= maxcutoff) or (gnomeAF_nfe >= mincutoff and  gnomeAF_nfe <= maxcutoff):
                        return True
                else:
                        return False


def main():
	if len(sys.argv) < 6:
		print("Usage: python preprocess_gzvcf.py filtered_input.vcf.gz famfile mincutoff maxcutoff outprefix")
		exit(1)
	
	infile = sys.argv[1]
	famfile = sys.argv[2]
	outprefix = sys.argv[5]
	mincutoff = float(sys.argv[3])
	maxcutoff = float(sys.argv[4])

	of = open(outprefix+"_mutations.tsv", "w")
	genes = []
	mut_rec = []

	famlines = [l.strip() for l in open(famfile)]
	cases = []
	controls = []
	for l in famlines:
		if l[0] != "#":
			x = l.split('\t')
			if x[5] == '1':
				controls.append(x[1])
			elif x[5] == '2':
				cases.append(x[1])

	of.write("#Gene\tmutid\tchg\tnchg\tppchg\trsid\ttop_csq\tcsq\timpact\tpolyphen\tsift\tcsq2\tgg_af\tgg_nfe_af\tge_af\tge_nfe_af\tac\tan\tPopMAF\tac_case\tan_case\tac_control\tan_control\tmutated_individuals\n")
	with gzip.open(infile,"rt") as fp:
		l = fp.readline().strip()
		counter = 0
		while(l):
			if l[0:2] == "##":
				l = fp.readline().strip()
				continue
			elif l[0] == '#':
				info_ind = -1
				indiv_start = -1
				header = l.split('\t')
				for i in range(len(header)):
					if "INFO" in header[i]:
						info_ind = i
					elif "FORMAT" in header[i]:
						indiv_start = i+1
				indivs = [v for v in header[indiv_start:]]
				N = len(indivs)
				print(N)
				l = fp.readline().strip()
			elif l[0] != '#':
				x = l.split('\t')
				#print(x)
				chr = x[0]
				pos = x[1]
				if chr+':'+pos not in mut_rec:
					mut_rec.append(chr+':'+pos)
				else:
					l = fp.readline().strip()
					continue
				rsid = x[2]
				ref = x[3]
				alt = x[4]
				print(chr, pos, ref, alt)
				gts = x[indiv_start:]
				hets = 0
				homs = 0
				missing = 0
				mutated_indiv = []
				ac_case = 0
				ac_control = 0
				m_case = 0
				m_control = 0
				for i in range(len(gts)):
					if gts[i] == '0/1' or gts[i] == '1/0':
						hets += 1
						mutated_indiv.append(indivs[i])
						if indivs[i] in cases:
							ac_case += 1
						elif indivs[i] in controls:
							ac_control += 1
					elif gts[i] == '1/1':
						homs += 1
						mutated_indiv.append(indivs[i])
						mutated_indiv.append(indivs[i])
						if indivs[i] in cases:
                                                        ac_case += 2
						elif indivs[i] in controls:
                                                        ac_control += 2
					elif gts[i] == './.':
						missing += 1
						if indivs[i] in cases:
                                                        m_case += 1
						elif indivs[i] in controls:
                                                        m_control += 1
				af = 0.0
				ac = 2*homs + hets
				an = 2*(N-missing)
				af = ac*1.0/an
				an_case = 2*(len(cases)-m_case)
				an_control = 2*(len(controls)-m_control)
				#if af > 0.5:
				#	af = 1-af
				mutid = chr+':'+pos
				mut = ref + '>' + alt
				if len(ref) == 1 and len(alt) == 1:
					chg = 'SNP'
				else:
					chg = 'INDEL'
				info = x[info_ind].split(';')
				#if '=' not in info:
				#	pass
				ann_ind = -1
				dbgnom_g = -1
				dbgnom_e = -1
				dbgnom_ng = -1
				dbgnom_ne = -1
				dbaaref = -1
				dbaaalt = -1
				dbaapos = -1
				gnom_g = -1
				gnom_ng = -1
				pphen_ind = -1
				sift_ind = -1
				gnomgAF = 0.0
				gnomgAF_nfe = 0.0
				gnomeAF = 0.0
				gnomeAF_nfe = 0.0
				if len(info) <= 1:
					counter += 1
					print(info, counter)
					l = fp.readline().strip()
					continue
				for i in range(len(info)):
					if "ANN=" in info[i]:
						ann_ind = i
						#print(ann_ind)
					elif "dbNSFP_gnomAD_genomes_POPMAX_AF" in info[i]:
						dbgnom_g = i
					elif "dbNSFP_gnomAD_exomes_POPMAX_AF" in info[i]:
						dbgnom_e = i
					elif "dbNSFP_gnomAD_genomes_NFE_AF" in info[i]:
						dbgnom_ng = i
					elif "dbNSFP_gnomAD_exomes_NFE_AF" in info[i]:
						dbgnom_ne = i
					elif "gnomAD_genome=" in info[i]:
						gnom_g = i
					elif "gnomAD_genome_nfe" in info[i]:
						gnom_ng = i
					elif "dbNSFP_Polyphen2_HVAR_pred" in info[i]:
						pphen_ind = i
					elif "dbNSFP_SIFT_pred" in info[i]:
						sift_ind = i
					elif "dbNSFP_aaref" in info[i]:
						dbaaref = i
					elif "dbNSFP_aaalt" in info[i]:
						dbaaalt = i
					elif "dbNSFP_aapos" in info[i]:
                                                dbaapos	= i
				#print(info[ann_ind])
				if pphen_ind != -1:
					pphen = info[pphen_ind].split('=')[1]
				if sift_ind != -1:
					sift = info[sift_ind].split('=')[1]
				ppchg = "."
				aaref = aaalt = aapos = ppchg = ann = '.'
				if dbaaref != -1:
					aaref = info[dbaaref].split('=')[1]
				if dbaaalt != -1:
					aaalt =	info[dbaaalt].split('=')[1]
				if dbaapos != -1:
					aapos =	info[dbaapos].split('=')[1]
				if aaref != '.':
					ppchg = aaref+aapos+aaalt
				if ann_ind != -1:
					ann = info[ann_ind].split('|')
					#print(ann)
					gene = ann[3]
					impact = ann[2]
					csq = ann[1]
					top_csq = ann[1].split('&')[0]
					print(top_csq)
					if top_csq not in incl:
						l = fp.readline().strip()
						continue
					csq2 = 'unknown'
					if top_csq in damaging or 'D' in pphen or 'P' in pphen or 'D' in sift:
						csq2 = 'damaging'
					elif 'B' in pphen and 'T' in sift:
						csq2 = 'neutral'
					elif top_csq == 'synonymous_variant':
						csq2 = 'synonymous'
				if dbgnom_g != -1:
					dbg = info[dbgnom_g].split('=')[1].split(',')[0]
				else:
					dbg = 0
				if dbgnom_e != -1:
					dbe = info[dbgnom_e].split('=')[1].split(',')[0]
				else:
					dbe = 0
				if dbgnom_ng != -1:
					dbng = info[dbgnom_ng].split('=')[1].split(',')[0]
				else:
					dbng = 0
				if dbgnom_ne != -1:
					dbne = info[dbgnom_ne].split('=')[1].split(',')[0]
				else:
					dbne = 0
				if gnom_g != -1:
					gg = info[gnom_g].split('=')[1].split(',')[0]
				else:
					gg = 0
				if gnom_ng != -1:
					gng = info[gnom_ng].split('=')[1].split(',')[0]
				else:
					gng = 0
				if dbe != '.' and dbe != 0:
					gnomeAF = float(dbe)
				if dbne != '.' and dbne != 0:
					gnomeAF_nfe = float(dbne)
				if dbg != '.' and dbg != 0:
					gnomgAF = float(dbg)
				else:
					if gg != '0':
						gnomgAF = float(gg)
				if dbng != '.' and dbng != 0:
					gnomgAF_nfe = float(dbng)
				else:
					if gng != '0':
						gnomgAF_nfe = float(gng)
				if len(mutated_indiv)>0:
					genes.append(gene)
					s = gene + '\t' + mutid + '\t' +  chg + '\t' + mut + '\t' + ppchg + '\t' + rsid + '\t' + top_csq + '\t' + csq + '\t' + impact + '\t' + pphen + '\t' + sift + '\t' + csq2 + '\t' + str(gnomgAF) + '\t' + str(gnomgAF_nfe) + '\t' + str(gnomeAF) + '\t' + str(gnomeAF_nfe) + '\t' + str(ac) + '\t' + str(an) + '\t' + str(af) + '\t' + str(ac_case) + '\t' + str(an_case) + '\t' + str(ac_control) + '\t' + str(an_control) +'\t'+';'.join(mutated_indiv) + '\n'
					of.write(s)
                    
			l = fp.readline().strip()
	
	of.close()
	
if __name__ == '__main__':
	main()

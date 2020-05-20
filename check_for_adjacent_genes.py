#!/usr/bin/python

""" Find whether fusion pair contains adjacent genes """

import pandas as pd

homo_gtf = pd.read_csv('../ensemble_hg38_r94/Homo_sapiens.GRCh38.94.gtf', sep='\t', skiprows=5, header=None) 
comm3 = pd.read_csv('fusions_common_to_2_programs.csv', sep='\t')
comm3.columns = ['up_gene','dw_gene']
geneA_loc = pd.DataFrame()
geneB_loc = pd.DataFrame()
pair_table = pd.DataFrame()

for i in comm3.index:
	geneA = comm3.iloc[i,0]
	geneB = comm3.iloc[i,1]

	df = homo_gtf[(homo_gtf[2] == 'gene') & (homo_gtf[8].str.contains(geneA))][[0,4,6]]
	df2 = homo_gtf[(homo_gtf[2] == 'gene') &(homo_gtf[8].str.contains(geneB))][[0,3,6]]
	
	df['GeneA_ID'] = geneA
	df2['GeneB_ID'] = geneB

	geneA_loc = pd.concat([geneA_loc, df])
	geneB_loc = pd.concat([geneB_loc, df2])
	
	
geneA_loc.columns = ['chrA','stopA','strandA','GeneA_ID']	
geneB_loc.columns = ['chrB','startB','strandB','GeneB_ID']
pair_table = pd.concat([geneA_loc.reset_index(drop=True), geneB_loc.reset_index(drop=True), comm3], axis=1)
	
pair_table_filt = pair_table[(pair_table['chrA'] == pair_table['chrB']) & (pair_table['strandA'] == pair_table['strandB'])].reset_index(drop=True)

for i in pair_table_filt.index:
	chrom = pair_table_filt.iloc[i,0]
	stop = pair_table_filt.iloc[i,1]
	start = pair_table_filt.iloc[i,5]
	if start-stop > 0:
		tmp = homo_gtf[(homo_gtf[0] == chrom) & (homo_gtf[2] == 'gene') & (homo_gtf[3] > stop) & (homo_gtf[3] < start)]
		if tmp.empty:
			print('fusion pairs are adjacent genes!')
		else:
			print('{0:>2} {1:<20} {2:<10} {3:<20} {4:<10} {5:<30} {6:<10}'.format(chrom, 
							pair_table_filt.loc[i,'GeneA_ID'], stop, 
							pair_table_filt.loc[i,'GeneB_ID'], start, 
							'no. of genes in between', len(tmp)),
							sep='\t', file=open('genes_bw_pairs.txt','a'))
	else:
		tmp = homo_gtf[(homo_gtf[0] == chrom) & (homo_gtf[2] == 'gene') & (homo_gtf[3] > start) & (homo_gtf[3] < stop)]
		if tmp.empty:
			print('fusion pairs are adjacent genes!')
		else:
			
			print('{0:>2} {1:<20} {2:<10} {3:<20} {4:<10} {5:<30} {6:<10}'.format(chrom, 
							pair_table_filt.loc[i,'GeneA_ID'], stop, 
							pair_table_filt.loc[i,'GeneB_ID'], start, 
							'no. of genes in between', len(tmp)),
							sep='\t', file=open('genes_bw_pairs.txt','a'))
				


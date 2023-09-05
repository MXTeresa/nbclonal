import argparse
import pandas as pd



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat')
	args = parser.parse_args()
	#args.dat = 'data/no_cut_trans_freq.csv'
	dat = pd.read_csv(args.dat)
	all_a = []
	for t  in [0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.0025]:
		all_a.append(dat.query('ref != var & valid & snv_qualified1 & snv_qualified2').reset_index(drop=True).assign(
					clonal = lambda k: ((k['freq1'] >= 1-t) & (k['freq2'] <= t)) | 
						((k['freq1'] <= t) & (k['freq2'] >= 1-t))).\
				groupby(['SPECID1', 'SPECID2']).\
				agg({'clonal': sum}).reset_index().assign(t=t))

		
	all_a = pd.concat(all_a)
	all_a = all_a.pivot(index=['SPECID1', 'SPECID2'], columns='t', values='clonal').reset_index().\
		rename({'SPECID1': 'ID1', 'SPECID2':'ID2'})
	all_a.to_csv('data/iav_data_processed.csv')


if __name__ == "__main__":
	run()


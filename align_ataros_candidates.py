import argparse 
import csv
import jiwer
import os
import sys
from utils import *

def clean_text(text_series):
	regex0 = re.compile(r"'(?=[A-z]+\s?{[^}]+})")
	regex1 = re.compile(r'{[^}]+}|-|_|\.\.|\||\*')
	regex2 = re.compile(r'(\(([0-9]x|\?\?)\))|\[[^\]]+\]|-|\|')
	regex3 = re.compile(r'\(\s*([^?0-9)][^)]+)\)')
	clean = text_series.apply(lambda x: re.sub(regex0,'',x))
	clean = clean.apply(lambda x: re.sub(regex1,'',x))
	clean = clean.apply(lambda x: re.sub(regex2,'',x))
	clean = clean.apply(lambda x: re.sub(regex3,r'\1',x))
	clean = [re.sub(r'\s+',' ',c.strip()) for c in clean]
	clean = [c for c in clean if c]
	return clean

def clean_df(x):
	if args.keep_bracketed:
		regex = re.compile(r'({[^}]+}|-|\||\*)|(\([^\)]+\)|-|\|)')
	else:
		regex = re.compile(r'({[^}]+}|-|\||\*)|(\([^\)]+\)|\[[^\]]+\]|-|\|)')
	return re.sub(r'\s+',' ',re.sub(regex,'',x).strip())

def print_as_table(out,
				   candidates,
				   exclude_deletions=False,
				   exclude_insertions=False,
				   binary=False):
	ref = out.references[0]
	hyp = out.hypotheses[0]
	alignment = out.alignments[0]
	data = []
	deletion_indices = []
	insertion_indices = []
	for ax in alignment:
		if ax.type in ['equal','substitute']:
			if ax.type == 'equal':
				achar = 1 if binary else ''
			else:
				achar = 0 if binary else 'S'
			for ref_idx , hyp_idx in zip(range(ax.ref_start_idx,ax.ref_end_idx),
										 range(ax.hyp_start_idx,ax.hyp_end_idx)):
				rx = ref[ref_idx]
				hx = hyp[hyp_idx]
				cx = dict(candidates.iloc[hyp_idx])
				cx.update({'ref':rx,'hyp':hx,'code':achar})
				data.append(cx)

		elif ax.type == 'delete':
			if exclude_deletions:
				continue
			hx = '*'
			achar = 0 if binary else 'D'
			for ref_idx in range(ax.ref_start_idx,ax.ref_end_idx):
				deletion_indices.append(ref_idx)
				rx = ref[ref_idx]
				cx = {'text':'','prob':'','xmin':'','xmax':''}
				cx.update({'ref':rx,'hyp':hx,'code':achar})
				data.append(cx)
			
		elif ax.type == 'insert':
			if exclude_insertions:
				continue
			rx = '*'
			achar = 0 if binary else 'I' 
			for hyp_idx in range(ax.hyp_start_idx,ax.hyp_end_idx):
				hx = hyp[hyp_idx]
				cx = dict(candidates.iloc[hyp_idx])
				cx.update({'ref':rx,'hyp':hx,'code':achar})
				data.append(cx)
				insertion_indices.append(hyp_idx)
		else:
			import pdb;pdb.set_trace()
	return pd.DataFrame(data), deletion_indices	, insertion_indices


def main():
	candidates_path = args.indir
	textgrid_path = args.tg_path
	for sheet in glob.glob(os.path.join(candidates_path,"*score.txt")):
		speakers = os.path.basename(sheet).split('-')[:3]
		try:
			order = int(speakers[-1].split('_')[1].split('.')[0])-1
		except Exception:
			raise ValueError('Malformed file name (should include 2 speakers).')
		speakers[-1] = speakers[-1].split('_')[0]
		textgrid_match = glob.glob(os.path.join(textgrid_path,'*'+'*'.join(speakers)+'*'))
		assert(len(textgrid_match)==1)
		dfs,speakers2 = load_textgrid_intervals(textgrid_match[0],target_tier='transcription')
		df = dfs[speakers2.index(speakers[order])]
		text = clean_text(df.text)
		try:
			candidates = pd.read_csv(sheet,sep='\t',header=None,quoting=csv.QUOTE_NONE)
			candidates.columns = ['text','prob','xmin','xmax']
		except Exception:
			raise ValueError(f'Error loading transcript: {sheet}')
		candidates = candidates[candidates.text.apply(clean_df) != '']
		candidates['text'] = candidates.text.apply(clean_df)
		candidates.reset_index(inplace=True,drop=True)
		hyp = ' '.join(candidates.text.values)#.lower()

		ref = ' '.join(text)
		ref = re.sub(r"(.)(?=('t)|('ll)|('re)|('m)|('s)|('d)|('ve))",r'\1 ',ref)
		ref = re.sub(r'flagpole','flag pole', ref)
		ref = re.sub(r'bed(?=bug)','bed ', ref)
		ref = re.sub(r"(.)(?=(\.|,|!|\?))",r'\1 ',ref)
		ref = re.sub(r'\s+',' ',ref).strip()
		out= jiwer.process_words([ref],[hyp])
		vis = jiwer.visualize_alignment(out).split('\n')
		table, deletions ,insertions = print_as_table(out,candidates,binary=True)

		if args.keep_bracketed:
			table.loc[table.text.str.startswith('[_'),'code'] =1

		# # try and combine subtokens for 
		table['subtokens edited'] = np.nan
		i = 0
		while i < table.shape[0]:
			row = table.iloc[i]
			if (row.ref != row.hyp) and row.ref.lower().endswith(row.hyp.lower()) and (len(row.hyp) > 0):
				# walk back	
				for j in range(1,6):
					cand = table.iloc[i-j:i+1]
					if ''.join(cand.hyp).lower() == row.ref.lower():
						print('here')
						# add in row and change
						lens = [len(c) for c in cand.hyp]
						# if there are empty spaces in ref just change them
						precedents = table.iloc[i-j:i].ref.tolist()
						count_stars = precedents.count('*')

						# there's some stars, so add one less rows than the amount of stars
						unstarred = pd.DataFrame([
							['',np.nan,np.nan,np.nan,word,'*',0] for word in [
									w for w in precedents if w!='*']],
									columns=['text', 'prob', 'xmin', 'xmax', 'ref', 'hyp', 'code'],
									index=list(range(i-j,i-count_stars)))
						split_text = [row.ref[sum(lens[:k]):sum(lens[:k])+lens[k]] for k in range(len(lens))]
						table.loc[i-j:i,'ref'] = split_text
						table.loc[i-j:i,'subtokens edited'] = 1
						for k in range(i-j,i+1):
							row = table.iloc[k]
							if row.ref == row.hyp:
								table.at[k,'code'] = 1
							else:
								table.at[k,'code'] = 0

						if unstarred.shape[0] > 0:
							table = pd.concat([table.iloc[:i-j], unstarred, table.iloc[i-j:]]).reset_index(drop=True)
							i +=  j - count_stars + 1						
						break
				# if you finish the loop then nothing happened
				i += 1
			else:
				i+=1
		# if you get through the loop you've either added rows or failed to match
		
		table[['ref','text','code','prob','xmin','xmax','subtokens edited']].to_csv(os.path.basename(sheet).split('.')[0]+'err.tsv',sep='\t',index=False)
		print(sheet,out.wer,sep='\n')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--indir',type=str,
			help='folder where Whisper outputs are stored')
	parser.add_argument('--keep_bracketed',action='store_true',
			help='whether to retain bracketed outputs, such as [_TT_360_] during alignment')
	parser.add_argument('--tg_path',type=str,
			help='folder where human transcript TextGrids are stored')
	args = parser.parse_args()
	main()



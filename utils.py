import glob
import math
import numpy as np
import pandas as pd
import re
from praatio import textgrid
import os
import subprocess
import parselmouth
from parselmouth.praat import call

'''general utility functions for interacting with ATAROS textgrids'''


def load_textgrid_intervals(textgrid_path,target_tier='word'):
	# Load the TextGrid file
	with open(textgrid_path,'r') as f:
		text = f.read().strip()
	splittext = re.split(r'"IntervalTier"\n',text)
	preamble = splittext[0]
	good_tiers = [s for s in splittext if s[:len(target_tier)+2]==f'"{target_tier}"']
	good_tiers = [re.sub(r'[^\n]+\n[^\n]+\n"sp"\n','',good_tier) for good_tier in good_tiers]
	try:
		text = '"IntervalTier"\n'.join([preamble,good_tiers[0],good_tiers[1]])+'\n'
	except Exception:
		import pdb;pdb.set_trace()
	with open('temp.TextGrid','w') as f:
		f.write(text)
	try:
		tg = textgrid.openTextgrid('temp.TextGrid',
							   False,duplicateNamesMode='rename')
	except Exception:
		import pdb;pdb.set_trace()
	speakers = re.findall(r'NW[FM][0-9]+',textgrid_path)
	tier_data = [[],[]]
	for tier in tg.tiers:
		if tier.name.startswith(target_tier):
			speaker = speakers[-1] if tier.name.endswith('_2') else speakers[0]
			speaker_idx = speakers.index(speaker)
			for interval in tier.entries:
				tier_data[speaker_idx].append({
					'start_time': interval.start, 
					'end_time': interval.end, 
					'text': interval.label,
					'speaker': speaker})
	# Create a DataFrame from the extracted data
	dfs = [pd.DataFrame(t) for t in tier_data]
	return dfs,speakers

def load_stance_intervals(stance_paths,target_tier='coarse',split=False): 
	data = [[],[]] 
	for i,path in enumerate(stance_paths): 
		speaker = (os.path.basename(path)).split('-')[0]
		# rewrite the textgrid
		with open(path,'r') as f:
			text = f.read().strip()
		splittext = re.split(r'\s*item \[[0-9]\]:\s*\n',text)
		preamble = splittext[0]
		good_tier = [s for s in splittext if re.match(r'\s*class = "IntervalTier"\s*\n\s*name = "'+target_tier+'"',s)][0]
		good_tier = re.sub(r'[^\n]+\n[^\n]+\n[^\n]+\n\s*text = "sp"\s*\n','',good_tier)
		with open('temp.TextGrid','w') as f:
			f.write(preamble+'\t\titem [1]:\n'+good_tier)
		tg = textgrid.openTextgrid('temp.TextGrid',False,duplicateNamesMode='rename')
		tier = tg.getTier('coarse')
		for interval in tier.entries:
			data[i].append({
					'start_time': interval.start,'end_time': interval.end, 
					'stance': interval.label,'speaker': speaker})
	if split:
		to_return1  = pd.DataFrame(data[0])
		to_return2 = pd.DataFrame(data[1])
		if to_return1.shape[0] <= 0:
			to_return1 = None
		if to_return2.shape[0] <= 0:
			to_return2 = None
		return to_return1, to_return2
	else:
		to_return = pd.DataFrame(data[0]+data[1])
		if to_return.shape[0] > 0:
			return to_return
		else:
			return None

def extract_wave_segment(input_file, output_file, start_time_sec, end_time_sec):
	command = [
		'sox',
		input_file,
		'temp.wav',
		'trim', str(start_time_sec), str(end_time_sec-start_time_sec)
	]
	# Run the command using subprocess
	subprocess.run(command, check=True)
	
	sound = parselmouth.Sound('temp.wav')
	sound.scale(0.99)
	sound.save(output_file, "WAV")


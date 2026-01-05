import os
import numpy as np


def csvread(fname, return_ax=True, return_dsc=True):
	# ------------------------------------------------------------------- #
	def parse_numeric(s: str):
		"""Return int or float if `s` is a numeric string; otherwise return the string."""
		try:
			if "." in s or "e" in s.lower() in s:
				return float(s)
			else:
				return int(s)
		except ValueError:
			return s

	if not os.path.isfile(fname):
		print(f"Error: {fname} does not exist or is not a file.")
		return None

	with open(fname, "r", encoding="utf-8") as f:  # , encoding="utf-8"
		raw_lines = f.read().splitlines()
	
	marker = None
	marker_idx = None
	MainDictionary = {}


	marker = ''
	dsc = {}
	ax = {'xlabel':'', 'ylabel':'', 
	   'freq1':0.0,
	   'mwpower':0.0,
	   'title':''}
	
	pointID = 0
	# newSection = False
	hasMeas = False

	for ii, raw in enumerate(raw_lines):
		line = raw.strip()
		# newSection = False
		if not line: continue

		if ';' not in line:        # it is a "section"
			marker = line
			if marker=='Meas':
				nPoints = len(raw_lines)-ii-2   # current line and next line are not part of the deal
				ax['x'] = np.zeros(nPoints, dtype=float)
				data = np.zeros(nPoints, dtype=float)
				pointID = 0
			continue
		else:
			parts = [p.strip() for p in line.split(";")]
			if marker=='Meas':
				hasMeas = True
				if parts[0][0].isdigit():
					ax['x'][pointID]=parse_numeric(parts[0])
					data[pointID]=parse_numeric(parts[1])
					pointID+=1
				else:
					ax['xlabel'] = parts[0]
					ax['ylabel'] = parts[1]
			else:
				key = parts[0]
				key1 = ''.join(filter(str.isalpha, key))
				value = parse_numeric(parts[1]) if len(parts) > 1 else None
				text = parts[2] if len(parts) > 2 else ''
				dsc[f'{marker}.{key1}']=[value,text]
	if not hasMeas:
		if return_ax and return_dsc:
			return None, None, None
		elif return_ax or return_dsc:
			return None, None
		else:
			return None

	if 'Additional.Frequency' in dsc.keys():
		ax['freq1']=dsc['Additional.Frequency'][0]
	if 'Recipe.MicrowavePower' in dsc.keys():
		ax['mwpower']=dsc['Recipe.MicrowavePower'][0]
	if '.Name' in dsc.keys():
		ax['title']=dsc['.Name'][0]
	if return_ax and return_dsc:
		return ax, data, dsc
	elif return_ax:
		return ax, data
	elif return_dsc:
		return data, dsc
	else:
		return data

if __name__ == "__main__":
	fname = 'd:\\PSU\\PSUDrive\\OneDrive - The Pennsylvania State University\\Papers\\EPR_CbHydA1\\figures\\20220322_160928183_CbHydA1_pH6_Hox_CO_40K_30dB_03_result.csv'
	ax, data, dsc = csvread(fname, return_ax=True, return_dsc=True)
	import matplotlib.pyplot as plt
	plt.plot(ax['x'], data, color='b')   # line plot
	plt.title(label=ax['title'])
	plt.xlabel(ax['xlabel'])
	plt.show()
	print("\n=== Parsed data ===")
	print(ax)
	print(dsc)

MODELS = {
	'$\Lambda$CDM':
	{
		'omegas': [0.3,0,0.7],
		'w': [-1],
		'line': ['solid'],
	},
	'$\Lambda$CDM, w = -0.85':
	{
		'omegas': [0.3,0,0.7],
		'w': [-0.85],
		'line': [(0, (1, 5)),'dotted'],
	},
	'$\Lambda$CDM, w = -1.15':
	{
		'omegas': [0.3,0,0.7],
		'w': [-1.15],
		'line': [(0, (1, 5)),'dotted'],
	},
	'E-deS':
	{
		'omegas': [1,0,0],
		'w': [0],
		'line': ['dashed'],
	},
	'Open Universe':
	{
		'omegas': [0.3,0.7,0],
		'w': [0],
		'line': ['dashdot','dashdot'],
	},
	'Q-$\Lambda$CDM':
	{
		'omegas': [0.3,0,0.7],
		'w': ['Linear-redshift'],#, 'Chevallier-Polarski-Linder'],#, 'Barboza-Alcaniz', 'Low Correlation', 'Jassal-Bagla-Padmanabhan'],
		'line': ['dashdot','dashdot'],
	},
	'$\Lambda$CDM, different w':
	{
		'omegas': [0.3,0,0.7],
		'w': [-0.85,-1.15],
		'line': [(0, (1, 5)),'dotted'],
	},
}


DATASETS = {
	'6dF': 
	{
		'name': "6dF_full.txt",
		'lines_to_skip': 0,
		'z_pos': 15,
		'convert': 'v_km',
		'degrees_area': 17000,
	},

	'V4':
	{
		'name': "VIPERS_W4_SPECTRO_PDR2.txt",
		'lines_to_skip': 1,
		'z_pos': 10,
		'convert': 'none',
		'degrees_area': 24,
	},

	'V1': 
	{
		'name': "VIPERS_W1_SPECTRO_PDR2.txt",
		'lines_to_skip': 1,
		'z_pos': 10,
		'convert': 'none',
		'degrees_area': 24,
	},

	'DEEP': 
	{
		'name': "DEEP.txt",
		'lines_to_skip': 0,
		'z_pos': 8,
		'convert': 'none',
		'degrees_area': 0.03,
	},

	'FORS':
	{
		'name': "FORS.tsv",
		'lines_to_skip': 34,
		'z_pos': 1,
		'convert': 'none',
		'degrees_area': 0.014,
	},
	
	'TKRS':
	{
		'name': "TKRS.tsv",
		'lines_to_skip': 33,
		'z_pos': 1,
		'convert': 'none',
		'degrees_area': 0.04,
	},

	'GDDS':
	{
		'name': "GDDS.tsv",
		'lines_to_skip': 33,
		'z_pos': 1,
		'convert': 'none',
		'degrees_area': 0.03,
	},

	'K20':
	{
		'name': "K20.tsv",
		'lines_to_skip': 33,
		'z_pos': 1,
		'convert': 'none',
		'degrees_area': 0.014,
	},
}


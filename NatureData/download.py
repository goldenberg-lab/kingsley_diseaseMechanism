import urllib.request

BASE_URL = 'http://datahub.stjude.org/datasets/MBSE2016/hg19/'

with open('trackDb.txt') as fp:
	url = ''
	file_name = ''
	for line in fp:

		if line.startswith('bigDataUrl'):
			print(line)
			url = line.split()[-1]
			url = BASE_URL + url
			file_name = url.split('/')[-1]

		if line.startswith('subGroups'):
			print(line)
			subgroup = (line.split()[1]).split('=')[1]
			print(subgroup)

			urllib.request.urlretrieve(url, subgroup+"_"+file_name)

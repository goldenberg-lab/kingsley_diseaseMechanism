'''
To parse and reorganize the enhancers from Nature data (data/nature16546-s2.csv)
Change 'chr1:10000-20000' to 'chr1	10000	20000' with order.
'''

def main():
	with open('data/nature16546-s2.csv') as fp:
		result_dictionary = {}
		for line in fp:
			tmp = line.strip().split(':')
			chromosome = tmp[0]

			tmp2 = tmp[1].split('-')
			start = tmp2[0]
			end = tmp2[1]

			if chromosome not in result_dictionary:
				result_dictionary[chromosome] = []

			result_dictionary[chromosome].append((start, end))

	with open('data/sorted_enhancers.txt', 'w') as op:
		for chr_num in range(1, 25):
			if chr_num == 23:
				chr_num = 'X'
			if chr_num == 24:
				chr_num = 'Y'

			lists = result_dictionary['chr'+str(chr_num)]
			lists.sort(key=lambda x: int(x[0]))

			print ('Writing chr%s with %d enhancers' % (chr_num, len(lists)))

			for enhancers in lists:
				op.write('%s\t%s\t%s\n' % ('chr'+str(chr_num), enhancers[0], enhancers[1]))

if __name__ == '__main__':
	main()
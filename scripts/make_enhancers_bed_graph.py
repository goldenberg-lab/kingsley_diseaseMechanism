

if __name__ == '__main__':
	with open('data/sorted_enhancers.txt') as fp, open('data/enhancers.bedGraph', 'w') as op:
		for line in fp:
			line = line.strip()
			op.write(line+"\t1.0\n")
from collections import Counter
import math

def calculate_entropy(data):
    counts = Counter(data)
    total_numbers = len(data)
    probabilities = [count / total_numbers for count in counts.values()]
    entropy = -sum(p * math.log2(p) for p in probabilities)
    if entropy!=0:
        pm=1/len(counts.keys())
        entropyal=-(pm * math.log2(pm))*len(counts.keys())
        return round(entropy, 5)/round(entropyal, 5)
    else:
        return entropy


with open('py_entrin.txt', 'r') as infile, open('outputx.txt', 'w') as outfile:
	#infile.readline()  # skip first line
	for line in infile:
		line = line.strip()
		list=line.split('\t')
		filtered_list = [x for x in list if x != '.']
		result_entropy = calculate_entropy(filtered_list)
		outfile.write(line + '\t' + str(result_entropy)+ '\n')
input_file = 'seq.in'
output_file = 'output.fasta'


with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    i=0
    for line in infile:
        i=i+1
        columns = line.strip().split()
        id_col = columns[1].split("@")[0]+"NN"+str(i)
        description_col = ",".join(columns[1:6]) # 第六列作为描述
        sequence_col = columns[6]  # 第七列作为序列
        outfile.write(f'>{id_col} {description_col}seq{sequence_col}\n')


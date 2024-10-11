import re

# Function to calculate the reverse complement of a DNA sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Function to calculate primer properties (length, GC content, melting temperature)
def calculate_primer_properties(primer):
    length = len(primer)
    gc_content = (primer.count('G') + primer.count('C')) / length * 100
    melting_temp = 64.9 + 41 * (primer.count('G') + primer.count('C') - 16.4) / length
    return length, gc_content, melting_temp

# Function to score secondary structures
def score_secondary_structures(primer):
    # Check for self-annealing and hairpin formation
    if re.search(r'(....).*\1', primer):  # Simple check for repeats
        return 0
    if re.search(r'(.)(.)(.)(.)\4\3\2\1', primer):  # Check for hairpins
        return 0
    return 1

# Function to score 3' end specificity
def score_3_end_specificity(primer):
    # Check for runs of repetitive nucleotides at the 3' end
    if primer[-4:] in ["GGGG", "TTTT"]:
        return 0
    return 1

# Function to check if a primer meets the design principles
def is_valid_primer(primer):
    length, gc_content, melting_temp = calculate_primer_properties(primer)
    secondary_structure_score = score_secondary_structures(primer)
    end_specificity_score = score_3_end_specificity(primer)

    if not (18 <= length <= 25):
        return False
    if not (50 <= melting_temp <= 65):
        return False
    if not (40 <= gc_content <= 65):
        return False
    if secondary_structure_score == 0:
        return False
    if end_specificity_score == 0:
        return False

    return True

# Step 1: Gather Input Information
input_file_path = r"D:\Desing primer\Primers_Output.txt"
output_file_path = r"D:\Desing primer\Primer_Design_Output.txt"
additional_output_file_path = r"D:\Desing primer\Additional_Primer_Design_Output.txt"

# Read sequences from the input file
with open(input_file_path, "r") as file:
    data = file.read()

# Regular expression to parse the input file
matches = re.findall(r"Region: (.*?)\nChromosome: (\w+), Start: (\d+), Stop: (\d+)\nFetched Sequence \(Length: \d+\):\n([A-Za-z]+)\n", data)

# Function to design primers
def design_primers(sequence, region_name):
    sequence_length = len(sequence)
    mid_point = sequence_length // 2

    # Define the regions for primer design
    up_region = sequence[:mid_point - 100]
    down_region = sequence[mid_point + 100:]

    primers_set1 = []
    primers_set2 = []

    # Design forward primer for the up region (Set 1)
    for i in range(0, len(up_region) - 20, 1):
        forward_primer_set1 = up_region[i:i+20]
        if is_valid_primer(forward_primer_set1):
            primers_set1.append((region_name + "_up_Fwd_Set1", forward_primer_set1))
            break

    # Design reverse primer for the down region (Set 1)
    for i in range(0, len(down_region) - 20, 1):
        reverse_primer_set1 = down_region[i:i+20]
        if is_valid_primer(reverse_primer_set1):
            primers_set1.append((region_name + "_down_Rev_Set1", reverse_complement(reverse_primer_set1)))
            break

    # Design forward primer for the up region (Set 2)
    for i in range(len(up_region) - 20, len(up_region) - 40, -1):
        forward_primer_set2 = up_region[i:i+20]
        if is_valid_primer(forward_primer_set2) and forward_primer_set2 != forward_primer_set1:
            primers_set2.append((region_name + "_up_Fwd_Set2", forward_primer_set2))
            break

    # Design reverse primer for the down region (Set 2)
    for i in range(len(down_region) - 20, len(down_region) - 40, -1):
        reverse_primer_set2 = down_region[i:i+20]
        if is_valid_primer(reverse_primer_set2) and reverse_primer_set2 != reverse_primer_set1:
            primers_set2.append((region_name + "_down_Rev_Set2", reverse_complement(reverse_primer_set2)))
            break

    return primers_set1, primers_set2

# Process each parsed region and write outputs
with open(output_file_path, "w") as output_file, open(additional_output_file_path, "w") as additional_output_file:
    for match in matches:
        region, chromosome, start, stop, sequence = match
        start = int(start)
        stop = int(stop)

        # Design primers
        primers_set1, primers_set2 = design_primers(sequence, region)

        # Write output for Set 1
        output_file.write(f"Region: {region} Set 1\n")
        output_file.write(f"Chromosome: {chromosome}, Start: {start}, Stop: {stop}\n")
        additional_output_file.write(f"Region: {region}\n")
        additional_output_file.write(f"Chromosome: {chromosome}, Start: {start}, Stop: {stop}\n")
        additional_output_file.write(f"Fetched Sequence (Length: {len(sequence)}):\n{sequence}\n")
        
        for primer in primers_set1:
            region_name, primer_seq = primer
            length, gc_content, melting_temp = calculate_primer_properties(primer_seq)
            secondary_structure_score = score_secondary_structures(primer_seq)
            end_specificity_score = score_3_end_specificity(primer_seq)
            output_file.write(f">{region_name}\n{primer_seq}\n")
            output_file.write(f"Length: {length}, GC Content: {gc_content:.2f}%, Melting Temp: {melting_temp:.2f}°C\n")
            output_file.write(f"Secondary Structure Score: {secondary_structure_score}, 3' End Specificity Score: {end_specificity_score}\n")
            
            additional_output_file.write(f">{region_name}\n{primer_seq}\n")
            additional_output_file.write(f"Length: {length}, GC Content: {gc_content:.2f}%, Melting Temp: {melting_temp:.2f}°C\n")
            additional_output_file.write(f"Secondary Structure Score: {secondary_structure_score}, 3' End Specificity Score: {end_specificity_score}\n")
        
        output_file.write("\n")
        additional_output_file.write("\n")

        # Write output for Set 2
        output_file.write(f"Region: {region} Set 2\n")
        output_file.write(f"Chromosome: {chromosome}, Start: {start}, Stop: {stop}\n")
        additional_output_file.write(f"Region: {region}\n")
        additional_output_file.write(f"Chromosome: {chromosome}, Start: {start}, Stop: {stop}\n")
        additional_output_file.write(f"Fetched Sequence (Length: {len(sequence)}):\n{sequence}\n")
        
        for primer in primers_set2:
            region_name, primer_seq = primer
            length, gc_content, melting_temp = calculate_primer_properties(primer_seq)
            secondary_structure_score = score_secondary_structures(primer_seq)
            end_specificity_score = score_3_end_specificity(primer_seq)
            output_file.write(f">{region_name}\n{primer_seq}\n")
            output_file.write(f"Length: {length}, GC Content: {gc_content:.2f}%, Melting Temp: {melting_temp:.2f}°C\n")
            output_file.write(f"Secondary Structure Score: {secondary_structure_score}, 3' End Specificity Score: {end_specificity_score}\n")
            
            additional_output_file.write(f">{region_name}\n{primer_seq}\n")
            additional_output_file.write(f"Length: {length}, GC Content: {gc_content:.2f}%, Melting Temp: {melting_temp:.2f}°C\n")
            additional_output_file.write(f"Secondary Structure Score: {secondary_structure_score}, 3' End Specificity Score: {end_specificity_score}\n")
        
        output_file.write("\n")
        additional_output_file.write("\n")

print(f"Output saved to: {output_file_path}")
print(f"Additional output saved to: {additional_output_file_path}")

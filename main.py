# --------------fastaAligner--------------------

# This tool takes a fasta file with multiple entry sequences as input
# and outputs a report.txt file. The tool uses the Needleman-Wunsch
# Global Alignment Algorithm to align every seqeunce against each other
# in the fasta file. The report file contains the alignment score, the
# percent identity metric, and the alignment itself.

import os 
import math



#  HELPER FUNCTIONS FIRST


# This function checks two sequences and if they having matches characters,
# returns a 1, and a -1 if otherwise


def match_checker(seq_1, seq_2, i, j):
    if seq_1[j - 1] == seq_2[i - 1]:
        return 1
    else:
        return -1
    

# This function checks the northern element of the matrix from the current 
# position and returns its value


def check_N(nw_matrix, gap_penalty, i, j):
    return nw_matrix[i - 1][j] + gap_penalty


# This function checks the western element of the matrix from the current 
# position and returns its value



def check_W(nw_matrix, gap_penalty, i, j):
    return nw_matrix[i][j - 1] + gap_penalty



# This function checks the northwest element of the matrix from the current 
# position and returns its value plus a 1 if the two characters at the 
# current position are the same, and -1 if not (using the match_checker function)


def check_NW(nw_matrix, seq_1, seq_2, i, j):
    return nw_matrix[i-1][j-1] + match_checker(seq_1, seq_2, i, j)



# NEEDLEMAN-WUNSCH GLOBAL ALIGNER FUNCTION



 # This is the main alignment function that takes two sequences as input and returns
 # alignment information. 


def align(seq_1, seq_2):

    gap_penalty = -2    # FUNNY STORY HERE: So I was having so much trouble with an
                        # infinite loop running in my program and I ABSOLUTELY COULD NOT
                        # figure out why it was doing that. 45 minutes or so later:
                        # TURNS OUT I accidentally set the gap penalty to positive 2 instead 
                        # of negative 2. YIPPEE. Let's just say I had to
                        # run out to get some more painkillers the next day. 




    # STEP 1: Create matrices using sequence lengths for Needleman-Wunsch algorithm. The nw_matrix
    # contains alignment score information using Needleman-Wunsch scoring rubric. The direction_matrix 
    # contains information for when we traceback through the matrix to acquire the alignment string


    nw_matrix = [[0 for i in range(len(seq_1) + 1)] for j in range(len(seq_2) + 1)]

    direction_matrix = [[0 for i in range(len(seq_1) + 1)] for j in range(len(seq_2) + 1)]


    # Initialize gap values in matrix


    for i in range(len(seq_1) + 1):
        nw_matrix[0][i] = i * gap_penalty

    for j in range(len(seq_2) + 1):
        nw_matrix[j][0] = j * gap_penalty


    # STEP 2: Populate matrices going character by character in the two input sequences

    
    for i in range(1, len(seq_2) + 1):
        for j in range(1, len(seq_1) + 1):

            # Needleman-Wunsch algorithm

            NW_value = check_NW(nw_matrix, seq_1, seq_2, i, j)
            N_value = check_N(nw_matrix, gap_penalty, i, j)
            W_value = check_W(nw_matrix, gap_penalty, i, j)


            nw_matrix[i][j] = max(NW_value, N_value, W_value)

            # Populate direction matrix depending on grid score achieved in Needleman-Wunsch step

            if nw_matrix[i][j] == check_NW(nw_matrix, seq_1, seq_2, i, j):
                direction_matrix[i][j] = 'NW'
            elif nw_matrix[i][j] == check_N(nw_matrix, gap_penalty, i, j):
                direction_matrix[i][j] = 'N'
            elif nw_matrix[i][j] == check_W(nw_matrix, gap_penalty, i, j):
                direction_matrix[i][j] = 'W'


    # The alignment score is gathered from the last element of the score matrix

    alignment_score = nw_matrix[len(seq_2)][len(seq_1)]

    

    # STEP 3: This step reads the cardinal matrix, containing the N, NW, and W information,
    # and creates alignment strings so as to visualize the alignments 

    # Create new indices for the string while-loop

    di = len(seq_2)
    dj = len(seq_1)


    # Initalize alignment strings

    aligned_1 = ''
    aligned_2 = ''
    match_bars = ''

    # While loop starts at the last element of the matrix and moves in 
    # the direction of the Needleman-Wunsch algorithm

    while di > 0 and dj > 0:

        # Read the element and determine what value is there. Depending on 
        # the information the align_1 and align_2 strings are built accordingly

        if direction_matrix[di][dj] == 'NW':
            aligned_1 = seq_1[dj - 1] + aligned_1
            aligned_2 = seq_2[di - 1] + aligned_2

            # This step gathers information for the "match_bars" string, which will 
            # show which characters match each other

            if match_checker(seq_1, seq_2, di, dj) == 1:
                match_bars = '|' + match_bars
            else:
                match_bars = ' ' + match_bars    

            # Change indices accordingly

            di -= 1
            dj -= 1

        # Repeat for remaining cardinal values

        elif direction_matrix[di][dj] == 'W':
            aligned_1 = seq_1[dj - 1] + aligned_1
            aligned_2 = '-' + aligned_2
            match_bars = ' ' + match_bars
            dj -= 1

        # Repeat again

        elif direction_matrix[di][dj] == 'N':
            aligned_1 = seq_1[dj - 1] + aligned_1
            aligned_2 = seq_2[di - 1] + aligned_2
            match_bars = ' ' + match_bars
            di -= 1


    return alignment_score, aligned_1, aligned_2, match_bars    # SO MUCH LOVELY INFORMATION HERE





# This function counts the point mutations total and returns the percent 
# similarity between the two input alignments (with indels)


def percent_identity(align_1, align_2):

    counter = 0

    # Count point mutations between the two sequences

    for i in range(len(align_1)):
        if align_1[i] != align_2[i]:
            counter += 1

    result = (len(align_1) - counter) / len(align_1) * 100
    result = round(result, 1)

    return result



# This function formats the alignment so that it wraps the text after a 
# number of characters. The function does not return a value, but rather
# writes the strings directly into the output file included as a function
# argument

def align_wrap_writer(output_file, align_1, align_2, match_bars, max_width):
    for i in range(0, len(align_1), max_width):
        output_file.write(align_1[i : i + max_width])
        output_file.write("\n")
        output_file.write(match_bars[i : i + max_width])
        output_file.write("\n")
        output_file.write(align_2[i : i + max_width])
        output_file.write("\n"*2)
        






# MAIN PROGRAM
# This section reads the fasta file containing the sequence that are to be aligned
# and puts them into a dictionary. Then the sequences are aligned against each other
# and this will be outputted in a report.txt file containing the alignment score for each
# alignment, the total normalized alignment score, and each alignment in printed 


# Set directory and designate input fasta file name


print("\nStarting Program")

dir = os.getcwd() + "/"
print("Current working directory: ", dir)
user_input = input("Input file: ")
file_dir = dir + user_input

# Open input file in read mode

input_file = open(file_dir, "r")

# Initialize dictonary for fasta information

entries = {}

# Initialize ID and sequence variables

curr_seq = ""
curr_ID = ""

# Loop thorugh every line and extract sequence information into dictionary

line = input_file.readline().strip()

while line:
    curr_ID = line[1:]
    curr_seq = ""
    line = input_file.readline().strip()
    while line and not line.startswith(">"):
        curr_seq += line
        line = input_file.readline().strip()
    entries[curr_ID] = curr_seq

# Close input file

input_file.close()


# Create output file called "report.txt"

output_file_dir = dir + "report.txt"
output_file = open(output_file_dir, "w")

keys_list = list(entries.keys())

# This long loop is the main program that reads the fasta file and applies the 
# align() and percent_identity() functions to the information from the 
# {entries} dictionary

score_total = 0
counter_align = 1
length_keys_list = len(keys_list)

for i in range(len(keys_list)):
    for j in range(i + 1, len(keys_list)):
        if keys_list[i] != keys_list[j]:

            # Create ID for alignment. Use colon to separate deflines from each sequence being aligned

            alignment_ID = keys_list[i] + " : " + keys_list[j]

            # Use the align() function to align two sequences

            score, aligned_1, aligned_2, match_bars = align(entries[keys_list[i]], entries[keys_list[j]])
            score_total += score

            

            # Call percent_identity() function 

            percent_ident = percent_identity(aligned_1, aligned_2)
            percent_ident_string = "Percent Identity: " + str(percent_ident) + "% \n" 
            alignment_score_string = "Alignment Score: " + str(score)

            # Write information to output file

            output_file.write(alignment_ID)
            output_file.write("\n")
            output_file.write(alignment_score_string)
            output_file.write("\n")
            output_file.write(percent_ident_string)
            output_file.write("\n")

            # Call wrapper function to print alignment so that it wraps
            # the text after 70 characters per string

            align_wrap_writer(output_file, aligned_1, aligned_2, match_bars, 70)

            output_file.write("\n")
            output_file.write("------------------------------------------")
            output_file.write("\n")

            print("Alignment", counter_align, "out of", math.comb(length_keys_list, 2), "completed")
            counter_align += 1

# Close output file

output_file.close()

print("Results contained in report.txt file")

print("Done")
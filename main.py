from Bio import SeqIO
import sys
import logging
import argparse
import os
from collections import Counter
import pandas
import numpy
import matplotlib.pyplot as plt


def read_commandline():
    """
    Command line arguments

    :return: argparse argument object
    """
    parser = argparse.ArgumentParser(description=f"QSUB submissiong for ONT data mapping to monkeypox reference")

    parser.add_argument('--input_fasta', '-i', required=True, help='aligned multifasta file.')
    parser.add_argument('--output_file', '-o', required=True,
                        help='summary output csv file.')

    args = parser.parse_args()
    return args

def check_arguments(args):
    """
    Check that paths provided exist and create output folder if it doesn't exist already.

    :param args: output from arg parse
    :return: status
    """
    if not os.path.exists(args.input_fasta):
        logging.error(f'File does not appear to exists: {args.input_fasta}. Please check')
        return 1
    return 0

def parse_fasta_file(input_fasta):
    """
    Check if the fasta file is as expected, contains '>' and fasta sequence
    Extract fasta sequences from fasta file
    :param fasta_file: fasta file path
    :return: list of dictionaries [ {'sample_id': str(),  'seq': str(sequence}, ... , {'sample_id': str(),  'seq': str(sequence} ]
    """
    with open(input_fasta, "r") as handle:
        fasta_sequences = SeqIO.parse(handle, "fasta")
        if any(input_fasta) is False:
            logging.error(f'File {input_fasta} does not contain fasta sequences in the expected format. '
                          f'Please check file and try again.')
            sys.exit()

    logging.info(f"File {input_fasta} contains sequences in the expected format.")

    all_seq_dict_list = []
    fasta_sequences = SeqIO.parse(open(input_fasta), 'fasta')
    for fasta in fasta_sequences:
            name, seq = fasta.id, str(fasta.seq)
            seq_dict = {"seq_id": name, "Sequence": seq}
            all_seq_dict_list.append(seq_dict)
    return all_seq_dict_list


def check_fasta_alignments(all_seq_dict_list):
    '''
    Check that the fastas have been aligned, are they all the same length
    '''
    lengths_of_sequences = []
    for fasta_dict in all_seq_dict_list:
        lengths_of_sequences.append(len(fasta_dict['Sequence']))
        
    average_seq_len = float(sum(lengths_of_sequences) / len(lengths_of_sequences))

    
    first_seq_length = float(len(all_seq_dict_list[0]['Sequence']))

    if average_seq_len != first_seq_length:
        print(f'Average sequence length ({average_seq_len}) is not equal to the length of the first sequence ({first_seq_length}). Implying that alignments are not all equal length.')
        sys.exit()
    else:
        print(f'Average sequence length ({average_seq_len}) is the same as first sequence length ({first_seq_length})')
    
    return

def count_calls_per_position(all_seq_dict_list):
    """
    List of dictionaries contaiing sequences, cycle through each sequence and per position add base to appropriate list.
    Consolidate lists into single dataframe
    """
    # Count number of positions
    number_of_samples = len(all_seq_dict_list)
    print(f'Calculating the the number of unique calls per genome position, across all {number_of_samples} samples.')
    list_of_positon_dictionaries = []
    sequence_length = len(all_seq_dict_list[0]['Sequence'])
    # Slow loop below, probably could be made quicker, maybe with pandas?
    for pos in range(0, int(sequence_length)):
        list_of_base_calls = [seq_dic['Sequence'][pos].upper() for seq_dic in all_seq_dict_list]
        count_dic = Counter(list_of_base_calls)
        genome_position = int(pos)+1
        count_dic['Genome_position'] = genome_position
        list_of_positon_dictionaries.append(count_dic)

    position_df = pandas.DataFrame(list_of_positon_dictionaries)
    return position_df

def calculate_row_sums(df, cols_to_sum, sum_col_name):
    """
    :params: df = dataframe to sum rows, cols_to_sum = list of column names to sum values of, sum_col_name = str(name of column for sum values)
    :return: df with sum values added
    """
    df[sum_col_name] = df.loc[:, cols_to_sum].sum(axis=1)
    return df

def format_dataframe(position_df, output_file):
    """
    Take in list of dictionaries
    Order table
    Output as table
    """
    # Check column names for non-IUPAC characters
    column_names = position_df.columns.values.tolist()
    expected_column_names_list = ['Genome_position']
    iupac_list = ['A','C','G','T','U','R','Y','S','W','K','M','B','D','H','V','N','.','-','X']
    expected_column_names_list.extend(iupac_list)
    non_iupac_characters = list(set(column_names).difference(expected_column_names_list))
    # Test if items any items in column names list are not in expected column names list
    if len(non_iupac_characters) > 0: # if the number of unique items in list 2 is greater than 1
        print(f"Error: identifed unexpected non-IUPAC codes in multi-fasta alignment file. Please check sequences for the following characters:\n{', '.join(non_iupac_characters)}")

    # Calculate total number of samples per pos
    column_names.remove('Genome_position') # Remove position number from comparison
    position_df = calculate_row_sums(position_df, column_names, str('total_call_counts'))

    # Calculate the number of nucleotide calls
    if 'U' in column_names: # Account for RNA
        position_df = calculate_row_sums(position_df, ['A','T','G','C','U'], str('nt_call_counts'))
    else:
        position_df = calculate_row_sums(position_df, ['A','T','G','C'], str('nt_call_counts'))

    mixed_bases = []
    for mixed_base in ['R','Y','S','W','K','M','B','D','H','V']:
        if mixed_base in column_names:
            mixed_bases.append(mixed_base)

    position_df = calculate_row_sums(position_df, mixed_bases, str('mixed_call_counts'))

    if 'X' in column_names:
        position_df['ambig_count'] = position_df['N']+position_df['X']
    else:
        position_df['ambig_count'] = position_df['N']   


    # Calculate percentages
    position_df['nt_perc'] = 100/position_df['total_call_counts'] * position_df['nt_call_counts']
    if 'X' in column_names:
        position_df['ambig_perc'] = 100/position_df['total_call_counts'] * position_df['ambig_count']
    else:
        position_df['ambig_perc'] = 100/position_df['total_call_counts'] * position_df['N']

    position_df['gaps_perc'] = 100/position_df['total_call_counts'] * position_df['-']
    position_df['mixed_perc'] = 100/position_df['total_call_counts'] * position_df['mixed_call_counts']

    position_df.to_csv(output_file, index=False)
    return position_df

##TODO: Plot out into graph?
def plot_to_graph(position_df_plus, output_file):
    print(position_df_plus.head())
    x = position_df_plus['Genome_position'].tolist()
    y1 = position_df_plus['nt_call_counts'].tolist()
    y2 = position_df_plus['mixed_call_counts'].tolist()
    y3 = position_df_plus['-'].tolist()
    y4 = position_df_plus['ambig_count'].tolist()
    
    plt.style.use('seaborn-pastel')
    plt.stackplot(x,y1, y2, y3, y4, labels=['Nt calls','Mixed base calls','Gaps', 'Ambiguous Calls'])
    plt.legend(loc='upper left')
    print(output_file)
    output_file_name =  str(output_file+'.png')
    print(output_file_name)
    plt.savefig(output_file_name, dpi='figure', format='png', metadata=None,
        bbox_inches=None, pad_inches=0.05,
        facecolor='auto', edgecolor='auto',
        backend=None)


def main(args):
      # Create paths for qsub submission
    output_file = os.path.abspath(args.output_file)
    input_fasta = os.path.abspath(args.input_fasta)

    # all_seq_dict_list = parse_fasta_file(input_fasta)
    # check_fasta_alignments(all_seq_dict_list)
    # position_df = count_calls_per_position(all_seq_dict_list)
    position_df = pandas.read_csv('test4.csv')
    position_df_plus = format_dataframe(position_df, output_file)
    plot_to_graph(position_df_plus, output_file)


if __name__ == '__main__':
    args = read_commandline()
    check = check_arguments(args)
    if check == 1:
        sys.exit(logging.error("Arguments provided were not expected. Please check log."))
    sys.exit(main(args))
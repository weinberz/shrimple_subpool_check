import re
import numpy as np
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from utils.file_handling import read_background_map

# TODO: correctly handle the background reference
BACKGROUND_REFERENCE = 'ATGAGGGAAAATAACCAGTCCTCTACACTGGAATTCATCCTCCTGGGAGTTACTGGTCAGCAGGAACAGGAAGATTTCTTCTACATCCTCTTCTTGTTCATTTACCCCATCACATTGATTGGAAACCTGCTCATCGTCCTAGCCATTTGCTCTGATGTTCGCCTTCACAACCCCATGTATTTTCTCCTTGCCAACCTCTCCTTGGTTGACATCTTCTTCTCATCGGTAACCATCCCTAAGATGCTGGCCAACCATCTCTTGGGCAGCAAATCCATCTCTTTTGGGGGATGCCTAACGCAGATGTATTTCATGATAGCCTTGGGTAACACAGACAGCTATATTTTGGCTGCAATGGCATATGATCGAGCTGTGGCCATCAGCCGCCCACTTCACTACACAACAATTATGAGTCCACGGTCTTGTATCTGGCTTATTGCTGGGTCTTGGGTGATTGGAAATGCCAATGCCCTCCCCCACACTCTGCTCACAGCTAGTCTGTCCTTCTGTGGCAACCAGGAAGTGGCCAACTTCTACTGTGACATTACCCCCTTGCTGAAGTTATCCTGTTCTGACATCCACTTTCATGTGAAGATGATGTACCTAGGGGTTGGCATTTTCTCTGTGCCATTACTATGCATCATTGTCTCCTATATTCGAGTCTTCTCCACAGTCTTCCAGGTTCCTTCCACCAAGGGCGTGCTCAAGGCCTTCTCCACCTGTGGTTCCCACCTCACGGTTGTCTCTTTGTATTATGGTACAGTCATGGGCACGTATTTCCGCCCTTTGACCAATTATAGCCTAAAAGACGCAGTGATCACTGTAATGTACACGGCAGTGACCCCAATGTTAAATCCTTTCATCTACAGTCTGAGAAATCGGGACATGAAGGCTGCCCTGCGGAAACTCTTCAACAAGAGAATCTCCTCGGGGTAG'

def get_backgrounds(filename, receptor='OR1A1'):
    """
    Get the genetic backgrounds for the analysis.
    
    Returns:
    - List of genetic backgrounds.
    """
    # Define the genetic backgrounds
    backgrounds = read_background_map(filename, receptor=receptor)

    background_list = [background for background in backgrounds['background_info']]
    background_map = {row.background_info: row.background_sequence for _, row in backgrounds.iterrows()}
    
    return background_list, background_map

def filter_amplicons(df, sizes):
        filtered_data = []

        size_ranges = []
        for size in sizes:
                lower_bound = size * 0.95
                upper_bound = size * 1.05
                size_ranges.append((lower_bound, upper_bound))

        for size_range in size_ranges:
                lower_bound, upper_bound = size_range
                filtered_sample_data = df[(df['length'] >= lower_bound) & (df['length'] <= upper_bound)]
                filtered_data.append(filtered_sample_data)
    
        # Combine all filtered data
        filtered_data = pd.concat(filtered_data, ignore_index=True)
        return filtered_data

def clean_data(df, sample_names, sample_pools, amplicon_sizes):
    cleaned_data = []
    for sample_id, sample in enumerate(sample_names):
        sample_data = df[df['sample'] == sample]
        amplicons_in_sample = [amplicon_sizes[amp-1] for amp in sample_pools[sample_id]]
        sample_data = filter_amplicons(sample_data, amplicons_in_sample)
        
        cleaned_data.append(sample_data)

    # Combine all cleaned data
    cleaned_data = pd.concat(cleaned_data, ignore_index=True)
    return cleaned_data

def label_sequences(df, samples, labels, sizes):
    label_sizes = []
    for sample_labels in labels:
        sample_sizes = [sizes[label - 1] for label in sample_labels]
        label_sizes.append(sample_sizes)
    for index, sample in enumerate(samples):
        df_filter = df['sample'] == sample
        for label, size in zip(labels[index], label_sizes[index]):
            lower_bound = size * 0.95
            upper_bound = size * 1.05
            df.loc[df_filter & (df['length'] >= lower_bound) & (df['length'] <= upper_bound), 'label'] = label
    df.dropna(subset=['label'], inplace=True)
    return df

def extract_variant_barcode(seq, background_barcode, barcode_offset=17, barcode_length=10):
    if len(seq) < 100:
        idx = seq.find(background_barcode)
    else:
        idx = seq.find(background_barcode, len(background_barcode))

    if idx == -1 or idx > 100:
        return None, None
    
    var_start = idx + barcode_offset
    var_end = var_start + barcode_length
    barcode_position = idx
    if var_end > len(seq):
        return None, None
    return seq[var_start:var_end], barcode_position

def handle_mutation_check(result, reference, allowed_mismatches=0):
    windowed_seq = result['sequence'][100:]
    expected_fragment = result['expected_fragment']

    if pd.isna(expected_fragment):
        match = re.search(r'(\d+)-(\d+)', result['variant_info'])
        if match:
            expected_fragment = reference[3*int(match.group(1)):3*int(match.group(2))]
        else:
            result.update({'mutation_found': False, 'reason': 'No break in variant info'})
            return result
    
    alignment = pairwise2.align.localxs(expected_fragment, windowed_seq, -1, -1, one_alignment_only=True)[0]
    mismatches = len(expected_fragment) - alignment.score
    if mismatches <= allowed_mismatches:
        result['mutation_found'] = True
    else:
        result.update({'mutation_found': False, 'reason': 'Expected fragment mismatch'})
    
    result['observed_fragment'] = windowed_seq
    result['expected_fragment'] = expected_fragment
    result['mismatches'] = mismatches
    return result

def parse_mutation(mutation, variant_info):
    if pd.isna(mutation):
        return {'type': 'wildtype'} if 'wildtype' in variant_info.lower() else {'type': 'unknown'}
    if mutation.endswith('*'):
        return {'type': 'stop', 'position': mutation[:-1]}
    if mutation.startswith('ins'):
        return {'type': 'insertion', 'sequence': mutation[3:]}
    if mutation.startswith('del'):
        return {'type': 'deletion', 'info': mutation[3:]}
    return {'type': 'substitution', 'amino_acid': mutation}

def orient_and_extract_barcodes(df, background_barcode, barcode_offset=17, barcode_length=10):
    df['variant_barcode'] = None
    for sequence in df['sequence']:
        # Check if 'CTCTTTCTCT' exists within the first 100 bases
        alignment = pairwise2.align.localxs(sequence[:120],
                                            background_barcode, -1, -1, 
                                            one_alignment_only=True)[0]
        if alignment.score < 10:
            # check to see if the sequence is a reverse complement
            reverse_sequence = Seq(sequence).reverse_complement()
            alignment = pairwise2.align.localxs(reverse_sequence[:120],
                                                background_barcode, -1, -1, 
                                                one_alignment_only=True)[0]
            if alignment.score == 10:
                index_to_update = df[df['sequence'] == sequence].index
                sequence = reverse_sequence
                if not index_to_update.empty:
                    df.at[index_to_update[0], 'sequence'] = str(sequence)
            else:
                continue
        start_idx = alignment.start
        var_start = start_idx + barcode_offset
        var_end = var_start + barcode_length
    
        if var_end > len(sequence):
            continue
        index_to_update = df[df['sequence'] == sequence].index
        if not index_to_update.empty:
            df.at[index_to_update[0], 'variant_barcode'] = sequence[var_start:var_end]

    return df

def analyze_sequences(df_labeled, barcode_map, reference_dna, background_barcode):
    barcode_dict = barcode_map.set_index('variant_sequence').to_dict(orient='index')

    print(f"Total barcodes in map: {len(barcode_dict)}")
    print(f"Example barcodes: {list(barcode_dict.keys())[:3]}")

    df_labeled = orient_and_extract_barcodes(df_labeled, background_barcode)

    results = []
    for idx, row in df_labeled.iterrows():
        seq, sample, label, variant_barcode = row['sequence'], row['sample'], row['label'], row['variant_barcode']
        result = {
            'id': idx, 'sample': sample, 'label': label,
            'barcode_found': False, 'mutation_found': None,
            'sequence_length': len(seq), 'sequence': seq,
            'variant_barcode': variant_barcode 
        }

        if not variant_barcode:
            result['reason'] = 'WT barcode not found or invalid position'
            results.append(result)
            continue

        result.update({'barcode_found': True})

        if variant_barcode not in barcode_dict:
            result['reason'] = 'Variant barcode not in map'
            results.append(result)
            continue

        mut_info = barcode_dict[variant_barcode]
        codon, mutation, position, expected_fragment = mut_info['codon'], mut_info['mutation'], mut_info['position'], mut_info['expected_fragment']
        result.update({
            'variant_info': mut_info.get('variant_info', 'Unknown'),
            'subpool': mut_info.get('subpool', np.nan),
            'mutation': mutation,
            'position': position,
            'codon': codon,
            'expected_fragment': expected_fragment
        })

        parsed = parse_mutation(mutation, result['variant_info'])
        result['mutation_type'] = parsed['type']

        result = handle_mutation_check(result, reference_dna)
        
        results.append(result)
    
    results_df = pd.DataFrame(results)
    return results_df
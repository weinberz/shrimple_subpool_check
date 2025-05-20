import numpy as np
import pandas as pd
import re
from io import StringIO

def read_fastqs(filenames):
    read_info = []

    for filename in filenames:
        
        file = StringIO(filename.read().decode('utf-8'))
        sample = filename.name.split('/')[-1].split('.')[0]

        while True:
            header = file.readline()
            if not header:
                break  # EOF
            sequence = file.readline().strip()
            file.readline()  # plus line
            file.readline()  # quality line

            info = {'sample': sample,
                    'sequence': sequence,
                    'length' : len(sequence)}
            read_info.append(info)

    return pd.DataFrame(read_info)

def read_background_map(map_file, receptor='OR1A1'):
    """
    Compile the background map from the given files.
    
    Parameters:
    - map_file: Path to the barcode map file.
    - receptor: Receptor type to filter the barcodes.
    
    Returns:
    - DataFrame containing the compiled background map.
    """
    
    # Read in the barcode map
    background_map = pd.read_csv(map_file, sep='\t', header=None, names=['background_info', 'background_sequence'])
    # Drop rows not for our variant of interest
    background_map = background_map[background_map['background_info'].str.contains(receptor)]
    print(background_map['background_info'])
    background_map['background_info'] = [info.split('_')[3] for info in background_map['background_info']]

    return background_map

def read_barcode_map(map_file, identifier_file, background, receptor = 'OR1A1'):
    """
    Compile the barcode map from the given files.
    
    Parameters:
    - map_file: Path to the barcode map file.
    - identifier_file: Path to the oligo identifiers file.
    - background: Genetic background to filter the barcodes.
    
    Returns:
    - DataFrame containing the compiled barcode map.
    """
    
    # Read in the barcode map
    barcode_map = pd.read_csv(map_file, sep='\t', header=None, names=['variant_info', 'variant_sequence'])

    # Drop rows not for our variant of interest
    barcode_map = barcode_map[((barcode_map['variant_info'].str.contains(background)) |
                               (barcode_map['variant_info'].str.contains('wt'))) &
                            (barcode_map['variant_info'].str.contains(receptor))]

    # Extract position, mutation, and subpool
    barcode_map['position'] = barcode_map['variant_info'].str.extract(r'_([A-Z])(\d+)')[1]
    barcode_map['mutation'] = barcode_map['variant_info'].str.extract(r'_([A-Z])(\d+)([A-Z\*]+)$')[2]
    barcode_map['subpool'] = barcode_map['variant_info'].str.extract(r'_(\d+-\d+)_')[0]

    # Convert subpool to numerical ranking
    barcode_map['subpool'] = barcode_map['subpool'].str.split('-').str[0].astype(int).rank(method='dense').astype(int)

    def fill_missing_mutations(barcode_map):
        def infer_mutation(row):
            if pd.notna(row['mutation']):
                return row['mutation']
            vi = row['variant_info']
            if 'del' in vi:
                m = re.search(r'del\d+', vi)
                return m.group(0) if m else np.nan
            elif 'ins' in vi:
                m = re.search(r'ins[ACGT]+', vi)
                return m.group(0) if m else np.nan
            else:
                return np.nan

        barcode_map['mutation'] = barcode_map.apply(infer_mutation, axis=1)
        return barcode_map

    barcode_map = fill_missing_mutations(barcode_map)
    barcode_map['name'] = [variant_info[-1] for variant_info in barcode_map['variant_info'].str.split('_')]

    def pick_preferred(group):
         preferred = group[group['variant_info'].str.contains(background)]
         wt_present = group[group['variant_info'].str.contains('wt')]
         if not preferred.empty:
             return preferred.iloc[[0]]
         elif not wt_present.empty:
             return wt_present.iloc[[0]]
         else:
             return pd.DataFrame()

    barcode_map = barcode_map.groupby('name',group_keys = False).apply(pick_preferred)
    
    oligo_identifiers = pd.read_csv(identifier_file, sep='\t', header=0, names=['variant_info', 'codon', 'expected_fragment'])

    wt_barcodes = barcode_map[barcode_map['name']=='wildtype'].copy()
    barcode_map = pd.merge(barcode_map, oligo_identifiers, on='variant_info', how='inner')
    barcode_map = pd.concat([barcode_map, wt_barcodes])
    print('EXECUTING!!!!')
    print(barcode_map)

    return barcode_map
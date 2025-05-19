from utils.file_handling import read_barcode_map, read_fastqs
from utils import data_handling
from pathlib import Path
import streamlit as st

def project_root():
    return Path(__file__).resolve().parents[1]

AMPLICON_SIZES = [1148, 1052, 944, 845, 737, 629, 521, 413, 305]

@st.cache_data
def process_fastqs(file_info, genetic_background):
    """
    Process the uploaded FASTQ files and perform QC analysis.
    
    Parameters:
    - file_info: Dictionary containing file metadata.
    - genetic_background: Selected genetic background.
    
    Returns:
    - DataFrame containing the QC results.
    """
    print('New run starts here!')
    background, background_barcode = genetic_background
    filenames = [filename.split('.')[0] for filename in file_info.keys()]
    pools = [info['pools'] for info in file_info.values()]

    reads = read_fastqs([info['file'] for info in file_info.values()])
    print(len(reads))

    cleaned_reads = data_handling.clean_data(
        reads,
        filenames,
        pools,
        AMPLICON_SIZES
    )
    print(len(cleaned_reads))

    barcode_map = read_barcode_map(
        project_root() / "barcode_map_oligo.tsv",
        project_root() / "oligo_identifiers.tsv",
        background,
        receptor='OR1A1'
    )

    labeled_reads = data_handling.label_sequences(
        cleaned_reads,
        filenames,
        pools,
        AMPLICON_SIZES
    )
    print(len(labeled_reads))
    print(len(reads))

    analyzed_reads = data_handling.analyze_sequences(
        labeled_reads,
        barcode_map,
        data_handling.BACKGROUND_REFERENCE,
        background_barcode
    )
    
    return analyzed_reads
import streamlit as st
import pandas as pd
import numpy as np

from utils.data_handling import get_backgrounds
from utils.qc import process_fastqs
import altair as alt

st.set_page_config(page_title="Shrimple GG1 QC", layout="wide")
st.title("Shrimple GG1 QC")

# Sidebar for inputs
st.sidebar.header("Input Options")

backgrounds, background_map = get_backgrounds("barcode_map_gene.tsv")

# Genetic background selection
genetic_background = st.sidebar.selectbox(
    "Select a genetic background:",
    backgrounds
)

# File uploader
uploaded_files = st.sidebar.file_uploader(
    "Upload FASTQ files (one or more):",
    type=["fastq"],
    accept_multiple_files=True
)

# Dictionary to hold file metadata
file_info = {}

if uploaded_files:
    for uploaded_file in uploaded_files:
        st.subheader(f"File: {uploaded_file.name}")

        pools_str = st.sidebar.text_input(
            f"Pools for {uploaded_file.name} (e.g., 3,5,8):",
            key=uploaded_file.name
        )

        pools = [int(p.strip()) for p in pools_str.split(',') if p.strip().isdigit()]
        file_info[uploaded_file.name] = {
            "pools": pools,
            "file": uploaded_file
        }

# You can now use file_info to continue your downstream QC pipeline
# Add processing trigger
if st.button("Run QC Analysis"):
    with st.spinner("Running analysis... this may take a few minutes."):
        # Mock processing - replace with your actual script call
        
        results_df = process_fastqs(
            file_info,
            (genetic_background, background_map[genetic_background])
        )

        found_fraction = (results_df
                      .groupby('sample')
                      .agg({'barcode_found': lambda x: np.sum(x == True)/len(x),
                            'mutation_found': lambda x: np.sum(x == True)/len(x),})
                            .reset_index()
                      .rename(columns={'barcode_found': 'Background Barcodes Found',
                                       'mutation_found': 'Variant Barcodes Found'}))
        
        found_fraction_folded = found_fraction.melt(
            id_vars='sample',
            value_vars=['Background Barcodes Found', 'Variant Barcodes Found'],
            var_name='Category',
            value_name='Value'
        )

        # Reads per file
        reads_per_file = results_df.groupby('sample').agg({'sequence': 'count'}).reset_index()
        reads_per_file.rename(columns={'sequence': 'Reads'}, inplace=True)
        reads_per_file['Reads'] = reads_per_file['Reads'].astype(int)
        
        # Example outputs
        st.success("QC analysis complete!")

        st.subheader("Reads per File")
        st.dataframe(reads_per_file)
        st.subheader("Barcodes Found per Sample")
        chart = alt.Chart(found_fraction_folded).mark_bar().encode(
            x=alt.X('sample:N', title='Sample', axis=alt.Axis(labelAngle=45, labelLimit=1000)),
            y=alt.Y('Value:Q', title='Fraction Found'),
            color=alt.Color('Category:N'),
            xOffset='Category:N',
            tooltip=['sample:N', 'Category:N', 'Value:Q'],
        ).properties(
            width=300,
            height=400
        )

        st.altair_chart(chart, use_container_width=True)

        st.subheader("Reads per Subpool")
        
        subpool_counts = (results_df
                          .query('barcode_found == True')
                          .groupby('subpool')
                          .agg({'mutation_found': lambda x: np.sum(x==True)/len(x)})
                          .reset_index()
                          .rename(columns={'mutation_found': 'Fraction Found'}))
        
        reads_per_subpool = (subpool_counts
                             .melt(id_vars=['subpool'],
                                   value_vars=['Fraction Found'],
                                   var_name='Category',
                                   value_name='Value'))
        
        chart = alt.Chart(reads_per_subpool).mark_bar().encode(
            x=alt.X('subpool:N', title='Subpool', axis=alt.Axis(labelAngle=45, labelLimit=1000)),
            y=alt.Y('Value:Q', title='Fraction Found'),
            color=alt.Color('subpool:N'),
            tooltip=['subpool:N', 'Category:N', 'Value:Q'],
        ).properties(
            width=300,
            height=400
        )

        results_df['variant_barcode'] = results_df['variant_barcode'].astype(str)
        st.altair_chart(chart, use_container_width=True)
        st.dataframe(results_df)

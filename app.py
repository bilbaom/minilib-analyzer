import streamlit as st
import sys
import os
import tempfile
import pandas as pd
import plotly.express as px
from collections import Counter

# Ensure current directory is in sys.path for Posit Connect Cloud deployment
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import functions from the user's script
from extract_barcodes_mini import (
    build_pattern, process_file, is_nnk, find_stop_codons, gc_content, max_homopolymer_run
)

st.set_page_config(page_title="Barcode Extractor", layout="wide")
st.title("🧬 Sequencing FastQ Barcode Extractor")

st.markdown("""
Upload your sequencing library FASTQ file(s) to extract barcodes, count occurrences, and generate a summary CSV. 
""")

# Sidebar for parameters
with st.sidebar:
    st.header("Parameters")
    left_flank = st.text_input("Left Flank", value="ttgcagagctca")
    right_flank = st.text_input("Right Flank", value="aatacagctccc")
    
    use_fixed_length = st.checkbox("Use fixed barcode length?", value=False)
    if use_fixed_length:
        length = st.number_input("Barcode Length", min_value=1, value=21)
        length_tol = st.number_input("Length Tolerance", min_value=0, value=0)
        min_len = length - length_tol
        max_len = length + length_tol
    else:
        min_len = st.number_input("Minimum Length", min_value=1, value=1)
        max_len = st.number_input("Maximum Length", min_value=1, value=200)
        length = None
        
    mismatches = st.number_input("Allowed Mismatches (Flanks)", min_value=0, value=0)
    
    check_nnk_toggle = st.checkbox("Check NNK / Stop Codons?", value=False)
    check_nnk = st.number_input("Barcode Length for NNK Check", min_value=1, value=21) if check_nnk_toggle else None

uploaded_files = st.file_uploader("Upload FASTQ file(s)", type=["fastq", "fq", "gz"], accept_multiple_files=True)

if st.button("Run Analysis", type="primary"):
    if not uploaded_files:
        st.error("Please upload at least one FASTQ file.")
    else:
        with st.spinner("Analyzing reads..."):
            pattern = build_pattern(left_flank, right_flank, min_len, max_len, mismatches)
            
            counter = Counter()
            stats = {"total_reads": 0, "reads_with_barcode": 0}
            
            # Use temporary directory for uploaded files
            with tempfile.TemporaryDirectory() as tmpdir:
                for uploaded_file in uploaded_files:
                    tmp_path = os.path.join(tmpdir, uploaded_file.name)
                    with open(tmp_path, "wb") as f:
                        f.write(uploaded_file.getbuffer())
                    
                    process_file(tmp_path, pattern, counter, stats)
            
            # Build Results DataFrame
            results_data = []
            
            if check_nnk is not None:
                for bc, count in counter.most_common():
                    homo_len, homo_seq = max_homopolymer_run(bc)
                    if len(bc) == check_nnk:
                        nnk_flag = "Yes" if is_nnk(bc) else "No"
                        stops = find_stop_codons(bc)
                        stop_flag = "Yes" if stops else "No"
                        stop_seq = ",".join(stops) if stops else ""
                    else:
                        nnk_flag = "NA"
                        stop_flag = "NA"
                        stop_seq = "NA"
                    
                    results_data.append({
                        "barcode_seq": bc,
                        "length": len(bc),
                        "read_counts": count,
                        "GC_content": gc_content(bc),
                        "Max_homopolymer_run": homo_len,
                        "Max_homopolymer_seq": homo_seq,
                        "NNK_pattern": nnk_flag,
                        "Stop_codon": stop_flag,
                        "Stop_codon_seq": stop_seq
                    })
            else:
                 for bc, count in counter.most_common():
                    homo_len, homo_seq = max_homopolymer_run(bc)
                    results_data.append({
                        "barcode_seq": bc,
                        "length": len(bc),
                        "read_counts": count,
                        "GC_content": gc_content(bc),
                        "Max_homopolymer_run": homo_len,
                        "Max_homopolymer_seq": homo_seq
                    })
                    
            df = pd.DataFrame(results_data)
            
            # Display stats
            st.success("Analysis Complete!")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Reads Processed", stats['total_reads'])
            col2.metric("Reads with Barcode", stats['reads_with_barcode'])
            col3.metric("Unique Barcodes Found", len(counter))
            
            st.divider()
            
            if not df.empty:
                st.subheader("Results Data")
                st.dataframe(df)
                
                # Download button
                csv = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Results CSV",
                    data=csv,
                    file_name='barcode_counts.csv',
                    mime='text/csv',
                )
                
                st.divider()
                st.subheader("Data Visualizations")
                
                tab1, tab2, tab3 = st.tabs(["Count Distribution", "Size Distribution", "Stop Codons"])
                
                with tab1:
                    st.markdown("### Top 50 Barcodes by Read Count")
                    top_df = df.head(50)
                    fig_counts = px.bar(top_df, x='barcode_seq', y='read_counts', 
                                        labels={'barcode_seq': 'Barcode Sequence', 'read_counts': 'Read Counts'})
                    st.plotly_chart(fig_counts, use_container_width=True)
                    
                with tab2:
                    st.markdown("### Barcode Size Distribution")
                    fig_size = px.histogram(df, x='length', nbins=50, 
                                            labels={'length': 'Barcode Length', 'count': 'Number of Unique Barcodes'})
                    st.plotly_chart(fig_size, use_container_width=True)
                    
                with tab3:
                    if check_nnk is not None:
                        st.markdown("### Stop Codon Percentage")
                        # Filter to only checked lengths
                        checked_df = df[df['Stop_codon'] != 'NA']
                        if not checked_df.empty:
                            stop_counts = checked_df['Stop_codon'].value_counts().reset_index()
                            stop_counts.columns = ['Has Stop Codon', 'Count']
                            fig_stops = px.pie(stop_counts, values='Count', names='Has Stop Codon',
                                               color='Has Stop Codon', color_discrete_map={'Yes': 'red', 'No': 'green'})
                            st.plotly_chart(fig_stops, use_container_width=True)
                            
                            yes_count = checked_df[checked_df['Stop_codon'] == 'Yes'].shape[0]
                            st.info(f"Stop Codon Percentage (in valid length barcodes): **{(yes_count/len(checked_df)*100):.2f}%**")
                        else:
                            st.info("No barcodes matched the specified NNK check length.")
                    else:
                        st.info("Stop codon checking was not enabled. Enable it in the sidebar parameters to see this chart.")
            else:
                st.warning("No barcodes were found matching the criteria.")

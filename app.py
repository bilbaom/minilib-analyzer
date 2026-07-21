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
            
            # Save results in session state so interaction with widgets doesn't reset the page
            st.session_state["analysis_results"] = {
                "df": df,
                "stats": stats,
                "unique_barcodes": len(counter),
                "check_nnk": check_nnk
            }

st.divider()

# Parameter and Script Logic Guide
with st.expander("📖 Analysis Logic & Parameter Guide", expanded=("analysis_results" not in st.session_state)):
    st.markdown("""
    ### 🧬 How the Analysis Works
    1. **Dual-Strand Scanning**: For every read in your FASTQ file(s), the algorithm scans both the **forward sequence** and its **reverse complement**. This captures barcodes regardless of sequencing orientation.
    2. **Flank Pattern Matching**: Searches for sequences framed by `[LEFT FLANK] + [BARCODE] + [RIGHT FLANK]`.
    3. **Feature Extraction**:
       - Aggregates barcode counts to rank unique sequences.
       - Calculates **GC Content %** for each barcode.
       - Detects **Max Homopolymer Runs** (e.g. `AAAAA`) to highlight potential sequencing or synthesis artifacts.
       - Performs open reading frame (ORF) checks (NNK pattern and Stop codons) if enabled.

    ---

    ### ⚙️ Parameter Reference Guide
    * **Left / Right Flanks**: The constant DNA sequences flanking the variable barcode (default: `ttgcagagctca` and `aatacagctccc`).
    * **Barcode Length (Fixed vs. Range)**:
      * **Fixed Length**: Recommended when your library has a known length design (e.g., 21 bp for a 7-codon library). Prevents spurious short matches.
      * **Min / Max Length**: Allows extracting variable-length barcodes within a custom range [min, max].
    * **Allowed Mismatches**: Number of tolerated mismatches/edits per flank sequence (requires the `regex` library).
    * **NNK & Stop Codon Check**:
      * **NNK Codons**: Checks if codons follow $N-N-K$ ($N = A,C,G,T$; $K = G,T$), commonly used in mutagenesis libraries to encode all 20 amino acids while eliminating 2 of 3 stop codons.
      * **Stop Codons**: Identifies in-frame stop codons ($TAA$, $TAG$, $TGA$).
    """)

# Display results if available in session_state
if "analysis_results" in st.session_state:
    res = st.session_state["analysis_results"]
    df = res["df"]
    stats = res["stats"]
    unique_barcodes = res["unique_barcodes"]
    check_nnk = res["check_nnk"]

    # Display stats
    st.success("Analysis Complete!")
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Reads Processed", stats['total_reads'])
    col2.metric("Reads with Barcode", stats['reads_with_barcode'])
    col3.metric("Unique Barcodes Found", unique_barcodes)
    
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
        
        tab1, tab2, tab3, tab4, tab5 = st.tabs(["Count Distribution", "Size Distribution", "GC Content", "Homopolymers", "NNK & Stop Codons"])
        
        with tab1:
            st.markdown("### Barcode Count Distribution (All Barcodes)")
            plot_df = df.copy().reset_index(drop=True)
            plot_df['Rank'] = range(1, len(plot_df) + 1)
            
            log_y = st.checkbox("Log scale for Read Counts", value=False, key="log_y_scale")
            
            fig_counts = px.line(
                plot_df,
                x='Rank',
                y='read_counts',
                hover_data=['barcode_seq', 'length', 'GC_content'],
                labels={'Rank': 'Barcode Rank (Sorted by Count)', 'read_counts': 'Read Count'},
                log_y=log_y,
                title='Read Count per Barcode (Rank-Abundance Curve)'
            )
            st.plotly_chart(fig_counts, use_container_width=True)
            
        with tab2:
            st.markdown("### Barcode Size Distribution (Weighted by Read Counts)")
            fig_size = px.histogram(
                df,
                x='length',
                y='read_counts',
                histfunc='sum',
                labels={'length': 'Barcode Length (bp)', 'read_counts': 'Total Read Count'},
                title='Length Distribution (Weighted by Read Counts)'
            )
            st.plotly_chart(fig_size, use_container_width=True)

        with tab3:
            st.markdown("### GC Content Distribution (Weighted by Read Counts)")
            fig_gc = px.histogram(
                df,
                x='GC_content',
                y='read_counts',
                histfunc='sum',
                nbins=50,
                labels={'GC_content': 'GC Content (%)', 'read_counts': 'Total Read Count'},
                title='GC Content Distribution (Weighted by Read Counts)'
            )
            st.plotly_chart(fig_gc, use_container_width=True)

        with tab4:
            st.markdown("### Most Common Homopolymers (> 2 bp)")
            # Filter homopolymers strictly larger than 2 bp
            homo_df = df[df['Max_homopolymer_run'] > 2].copy()
            if not homo_df.empty:
                total_unique_barcodes = len(df)
                total_reads_with_barcodes = df['read_counts'].sum()

                homo_summary = homo_df.groupby("Max_homopolymer_seq").agg(
                    Run_Length=("Max_homopolymer_run", "first"),
                    Unique_Barcodes=("barcode_seq", "count"),
                    Total_Read_Count=("read_counts", "sum")
                ).reset_index()
                
                homo_summary["Read_Count_Pct"] = (homo_summary["Total_Read_Count"] / total_reads_with_barcodes) * 100
                homo_summary["Unique_Barcodes_Pct"] = (homo_summary["Unique_Barcodes"] / total_unique_barcodes) * 100
                
                homo_summary = homo_summary.sort_values(by="Total_Read_Count", ascending=False).reset_index(drop=True)
                
                top_homo = homo_summary.head(20)
                fig_homo = px.bar(
                    top_homo,
                    x='Max_homopolymer_seq',
                    y='Read_Count_Pct',
                    hover_data=['Run_Length', 'Total_Read_Count', 'Unique_Barcodes'],
                    labels={
                        'Max_homopolymer_seq': 'Homopolymer Sequence',
                        'Read_Count_Pct': '% of Total Reads',
                        'Total_Read_Count': 'Total Read Count',
                        'Unique_Barcodes': 'Unique Barcodes Count'
                    },
                    title='Top Most Frequent Homopolymers > 2 bp (% of Total Read Counts)',
                    color='Run_Length',
                    color_continuous_scale='Viridis'
                )
                st.plotly_chart(fig_homo, use_container_width=True)
                
                st.markdown("#### Summary Table of Homopolymer Sequences (> 2 bp)")
                display_summary = homo_summary.copy()
                display_summary['Read_Count_Pct'] = display_summary['Read_Count_Pct'].map(lambda x: f"{x:.2f}%")
                display_summary['Unique_Barcodes_Pct'] = display_summary['Unique_Barcodes_Pct'].map(lambda x: f"{x:.2f}%")
                
                st.dataframe(display_summary, use_container_width=True)
            else:
                st.info("No homopolymers larger than 2 bp detected.")

        with tab5:
            if check_nnk is not None:
                checked_df = df[df['NNK_pattern'] != 'NA']
                if not checked_df.empty:
                    total_unique = len(checked_df)
                    total_reads = checked_df['read_counts'].sum()
                    
                    nnk_yes_unique = (checked_df['NNK_pattern'] == 'Yes').sum()
                    nnk_yes_reads = checked_df[checked_df['NNK_pattern'] == 'Yes']['read_counts'].sum()
                    nnk_pct_u = (nnk_yes_unique / total_unique) * 100
                    nnk_pct_r = (nnk_yes_reads / total_reads) * 100 if total_reads > 0 else 0
                    
                    stop_yes_unique = (checked_df['Stop_codon'] == 'Yes').sum()
                    stop_yes_reads = checked_df[checked_df['Stop_codon'] == 'Yes']['read_counts'].sum()
                    stop_pct_u = (stop_yes_unique / total_unique) * 100
                    stop_pct_r = (stop_yes_reads / total_reads) * 100 if total_reads > 0 else 0

                    st.markdown("### NNK Pattern & Stop Codon Analysis")
                    col_m1, col_m2 = st.columns(2)
                    with col_m1:
                        st.metric("NNK Pattern Match (Unique Barcodes)", f"{nnk_pct_u:.2f}%", f"{nnk_yes_unique} / {total_unique}")
                        st.metric("NNK Pattern Match (Weighted by Reads)", f"{nnk_pct_r:.2f}%", f"{nnk_yes_reads} / {total_reads} reads")
                    with col_m2:
                        st.metric("Stop Codon Presence (Unique Barcodes)", f"{stop_pct_u:.2f}%", f"{stop_yes_unique} / {total_unique}")
                        st.metric("Stop Codon Presence (Weighted by Reads)", f"{stop_pct_r:.2f}%", f"{stop_yes_reads} / {total_reads} reads")

                    st.divider()
                    col_p1, col_p2 = st.columns(2)
                    with col_p1:
                        st.markdown("#### NNK Pattern Conformance")
                        nnk_counts = checked_df['NNK_pattern'].value_counts().reset_index()
                        nnk_counts.columns = ['Conforms to NNK', 'Count']
                        fig_nnk = px.pie(
                            nnk_counts, values='Count', names='Conforms to NNK',
                            color='Conforms to NNK', color_discrete_map={'Yes': 'green', 'No': 'red'}
                        )
                        st.plotly_chart(fig_nnk, use_container_width=True)

                    with col_p2:
                        st.markdown("#### Stop Codon Distribution")
                        stop_counts = checked_df['Stop_codon'].value_counts().reset_index()
                        stop_counts.columns = ['Has Stop Codon', 'Count']
                        fig_stops = px.pie(
                            stop_counts, values='Count', names='Has Stop Codon',
                            color='Has Stop Codon', color_discrete_map={'Yes': 'red', 'No': 'green'}
                        )
                        st.plotly_chart(fig_stops, use_container_width=True)
                else:
                    st.info("No barcodes matched the specified NNK check length.")
            else:
                st.info("NNK & Stop codon checking was not enabled. Enable it in the sidebar parameters to see this analysis.")
    else:
        st.warning("No barcodes were found matching the criteria.")

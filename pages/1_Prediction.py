import os
import sys
from types import ModuleType

m = ModuleType("aaindex._aaindex_matrix")

class MockMap(dict):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class MockAAIndex:
    
    def __init__(self, name):
        self.aaindex_file = name
        self.aaindex_json = {}
    def parse_aaindex(self):
        return {}

m.Map = MockMap
m.AAIndex = MockAAIndex

sys.modules["aaindex._aaindex_matrix"] = m

import pandas as pd
import protpy
import re
import joblib
import plotly.express as px
from Bio import SeqIO
from io import StringIO

st.set_page_config( page_title="AMP-PPI", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h1 style='text-align:center;'> AMP-PPI: A Machine Learning based tool for prediction of Protein-Protein Interaction between Antimicrobial Peptides and their Target Proteins</h1>", unsafe_allow_html=True)


st.markdown("---")
st.text("")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
HELIX_PROPENSITY = {"A": 1.45, "C": 0.77, "D": 0.98, "E": 1.53, "F": 1.12, "G": 0.53, "H": 1.24, "I": 1.00, "K": 1.07, "L": 1.34, "M": 1.20, "N": 0.73, "P": 0.59, "Q": 1.17, "R": 0.79, "S": 0.79, "T": 0.82, "V": 1.14, "W": 1.14, "Y": 0.61}
SHEET_PROPENSITY = {"A": 0.97, "C": 1.30, "D": 0.80, "E": 0.26, "F": 1.28, "G": 0.81, "H": 0.71, "I": 1.60, "K": 0.74, "L": 1.22, "M": 1.67, "N": 0.65, "P": 0.62, "Q": 1.23, "R": 0.90, "S": 0.72, "T": 1.20, "V": 1.65, "W": 1.19, "Y": 1.29}
COIL_PROPENSITY = {"A": 0.72, "C": 1.17, "D": 1.41, "E": 0.80, "F": 0.58, "G": 1.64, "H": 1.00, "I": 0.51, "K": 1.01, "L": 0.57, "M": 0.67, "N": 1.68, "P": 1.91, "Q": 0.98, "R": 1.24, "S": 1.33, "T": 1.03, "V": 0.50, "W": 0.99, "Y": 1.29}

def clean_sequence(seq):
    
    seq = str(seq).upper().strip().replace(" ", "")
    seq = re.sub(r'[^A-Z]', '', seq)
    seq = ''.join([aa for aa in seq if aa in VALID_AA])
    return seq if len(seq) > 0 else None

def extract_features(sequence: str):
    
    def to_dict_safe(x):
        if isinstance(x, pd.Series): return x.to_dict()
        elif isinstance(x, pd.DataFrame): return x.iloc[0].to_dict()
        return x

    n = len(sequence)
    ssp = {
        "Helix": round(sum(HELIX_PROPENSITY[aa] for aa in sequence) / n, 3),
        "Sheet": round(sum(SHEET_PROPENSITY[aa] for aa in sequence) / n, 3),
        "Coil": round(sum(COIL_PROPENSITY[aa] for aa in sequence) / n, 3)
    }

    row = {}
    row.update(to_dict_safe(protpy.amino_acid_composition(sequence)))
    row.update(to_dict_safe(protpy.conjoint_triad(sequence)))
    row.update(to_dict_safe(protpy.ctd_composition(sequence)))
    row.update(to_dict_safe(protpy.ctd_transition(sequence)))
    row.update(to_dict_safe(protpy.ctd_distribution(sequence)))
    row.update(ssp)
    return row

def run_prediction(combined_df):
    
    if combined_df.empty:
        st.warning("No valid sequences found within the specified length ranges.")
        return None
    try:
        model = joblib.load("static/model/svm_model.joblib")
        X = combined_df.drop(columns=['amp', 'protein'], errors='ignore')
        
        preds = model.predict(X)
        probs = model.predict_proba(X)[:, 1]
        
        results = combined_df[['amp', 'protein']].copy()
        results['Interaction_Score'] = probs
        results['Prediction'] = ["Positive" if p == 1 else "Negative" for p in preds]
        return results
    except Exception as e:
        st.error(f"Prediction Error: {e}")
        return None

def process_fasta_buffer(buffer, min_l, max_l):
    
    records = []
    for r in SeqIO.parse(buffer, "fasta"):
        cleaned = clean_sequence(r.seq)
        if cleaned and min_l <= len(cleaned) <= max_l:
            records.append({"id": r.id, "seq": cleaned})
    return records

st.title("🧬 AMP-Protein Interaction Predictor")

st.info("💡 New here? Download example files to see how the data should be formatted.")
c_dl1, c_dl2 = st.columns(2)
with c_dl1:
    if os.path.exists("static/data/example_amp.fasta"):
        with open("static/data/example_amp.fasta", "rb") as f:
            st.download_button("📂 Download Example AMPs (FASTA)", f, "example_amps.fasta")
with c_dl2:
    if os.path.exists("static/data/example_protein.fasta"):
        with open("static/data/example_protein.fasta", "rb") as f:
            st.download_button("📂 Download Example Proteins (FASTA)", f, "example_proteins.fasta")

st.markdown("---")

st.subheader("🛠️ Preprocessing Filters")
col_s1, col_s2 = st.columns(2)
with col_s1:
    amp_range = st.slider("AMP Length Filter", 5, 200, (5, 100), help="Strict length check for AMPs.")
with col_s2:
    prot_range = st.slider("Protein Length Filter", 50, 5000, (50, 2000), help="Strict length check for Target Proteins.")

st.markdown("---")

tab1, tab2, tab3 = st.tabs(["⌨️ Single Pair", "📂 Bulk Combinations", "🔗 Paired Upload"])

with tab1:
    st.subheader("Manual Input")
    c1, c2 = st.columns(2)
    with c1:
        s_amp = st.text_area("AMP Sequence", placeholder="KLSPSL...")
    with c2:
        s_prot = st.text_area("Protein Sequence", placeholder="MRPPQC...")
    
    if st.button("Predict Single Interaction", type="primary"):
        c_amp = clean_sequence(s_amp)
        c_prot = clean_sequence(s_prot)
        
        if c_amp and c_prot:
            
            if (amp_range[0] <= len(c_amp) <= amp_range[1]) and (prot_range[0] <= len(c_prot) <= prot_range[1]):
                f_amp = extract_features(c_amp)
                f_prot = extract_features(c_prot)
                comb_row = {"amp": "Manual_AMP", "protein": "Manual_Protein"}
                for k, v in f_amp.items(): comb_row[f"{k}_x"] = v
                for k, v in f_prot.items(): comb_row[f"{k}_y"] = v
                
                res = run_prediction(pd.DataFrame([comb_row]))
                if res is not None:
                    st.metric("Interaction Score", f"{res['Interaction_Score'].values[0]:.4f}")
                    st.success(f"Result: {res['Prediction'].values[0]}")
            else:
                st.error("Sequences do not meet the length requirements set in the sliders above.")

with tab2:
    st.subheader("Discovery Mode (All-vs-All)")
    u1, u2 = st.columns(2)
    with u1:
        file_amp = st.file_uploader("Upload amp.fasta", type=["fasta"], key="bulk_u_amp")
    with u2:
        file_prot = st.file_uploader("Upload protein.fasta", type=["fasta"], key="bulk_u_prot")
    
    use_example_bulk = st.checkbox("Use example sequences instead")
    
    if st.button("Run Bulk Analysis", type="primary"):
        if use_example_bulk or (file_amp and file_prot):
            with st.spinner("Processing..."):
                
                buf_a = open("static/data/example_amp.fasta", "r") if use_example_bulk else StringIO(file_amp.getvalue().decode())
                buf_p = open("static/data/example_protein.fasta", "r") if use_example_bulk else StringIO(file_prot.getvalue().decode())
                
                amps = process_fasta_buffer(buf_a, amp_range[0], amp_range[1])
                prots = process_fasta_buffer(buf_p, prot_range[0], prot_range[1])
                
                rows = []
                for a in amps:
                    f_a = extract_features(a['seq'])
                    for p in prots:
                        f_p = extract_features(p['seq'])
                        combined = {"amp": a['id'], "protein": p['id']}
                        combined.update({f"{k}_x": v for k, v in f_a.items()})
                        combined.update({f"{k}_y": v for k, v in f_p.items()})
                        rows.append(combined)
                
                st.session_state['results'] = run_prediction(pd.DataFrame(rows))

with tab3:
    st.subheader("🔗 Paired Prediction (1-to-1)")
    st.caption("Upload separate files to predict specific pairs (AMP 1 with Protein 1, etc.) instead of all combinations.")
    
    u3_col1, u3_col2 = st.columns(2)
    with u3_col1:
        file_amp_p = st.file_uploader("Upload amp.fasta", type=["fasta"], key="paired_u_amp")
    with u3_col2:
        file_prot_p = st.file_uploader("Upload protein.fasta", type=["fasta"], key="paired_u_prot")
    
    # Example trigger for this specific tab
    run_example_paired = st.button("🚀 Run Example Paired Prediction")

    if st.button("Predict Pairs", type="primary") or run_example_paired:
        # Determine source: Example files or User Uploads
        if run_example_paired:
            buf_a = open("static/data/example_amp.fasta", "r")
            buf_p = open("static/data/example_protein.fasta", "r")
        elif file_amp_p and file_prot_p:
            buf_a = StringIO(file_amp_p.getvalue().decode())
            buf_p = StringIO(file_prot_p.getvalue().decode())
        else:
            st.warning("Please upload both files or click 'Run Example'.")
            buf_a = buf_p = None

        if buf_a and buf_p:
            with st.spinner("Processing paired sequences..."):
                
                amps = process_fasta_buffer(buf_a, amp_range[0], amp_range[1])
                prots = process_fasta_buffer(buf_p, prot_range[0], prot_range[1])
                
                if run_example_paired:
                    buf_a.close()
                    buf_p.close()

                # Determine the number of pairs based on the shorter list
                num_pairs = min(len(amps), len(prots))
                
                if num_pairs == 0:
                    st.error("No valid pairs found within the specified length ranges.")
                else:
                    rows = []
                    for i in range(num_pairs):
                        
                        a_item = amps[i]
                        p_item = prots[i]
                        
                        f_a = extract_features(a_item['seq'])
                        f_p = extract_features(p_item['seq'])
                        
                        combined = {
                            "amp": f"{a_item['id']} ({a_item['seq'][:10]}...)",
                            "protein": f"{p_item['id']} ({p_item['seq'][:10]}...)"
                        }
                        
                        combined.update({f"{k}_x": v for k, v in f_a.items()})
                        combined.update({f"{k}_y": v for k, v in f_p.items()})
                        rows.append(combined)

                    st.session_state['results'] = run_prediction(pd.DataFrame(rows))
                    
                    if num_pairs < max(len(amps), len(prots)):
                        st.warning(f"File lengths differ. Only the first {num_pairs} pairs were processed.")

if 'results' in st.session_state and st.session_state['results'] is not None:
    df = st.session_state['results']
    st.markdown("---")
    st.subheader("📊 Results & Analytics")
    
    c_g1, c_g2 = st.columns([2, 1])
    with c_g1:
        fig = px.histogram(df, x="Interaction_Score", color="Prediction", 
                           title="Interaction Probability Distribution",
                           color_discrete_map={"Positive":"#00CC96", "Negative":"#EF553B"},
                           marginal="box")
        st.plotly_chart(fig, use_container_width=True)
    
    with c_g2:
        pie_fig = px.pie(df, names="Prediction", title="Summary of Interactions",
                         color="Prediction", color_discrete_map={"Positive":"#00CC96", "Negative":"#EF553B"})
        st.plotly_chart(pie_fig, use_container_width=True)

    st.dataframe(df.style.background_gradient(subset=['Interaction_Score'], cmap='RdYlGn'), use_container_width=True)
    
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Predictions CSV", csv, "results.csv", "text/csv")

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)

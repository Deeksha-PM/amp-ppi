import streamlit as st

st.set_page_config( page_title="AMP-PPI", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h1 style='text-align:center;'> AMP-PPI: A Machine Learning based tool for prediction of Protein-Protein Interaction between Antimicrobial Peptides and their Target Proteins</h1>", unsafe_allow_html=True)



st.markdown("---")
st.text("")

col1_1, col2_1 = st.columns([2, 2])

with col1_1:
    st.header("Background")
    st.text("Protein–protein interactions are central to understanding how antimicrobial peptides (AMPs) exert their effects on target proteins. Gaining insight into these interactions is important for exploring their applications in drug development and antimicrobial strategies. In this context, AMP-PPI, a machine learning-based predictor, has been developed to enable the prediction of protein–protein interactions (PPIs) between AMPs and their target proteins using amino acid sequence information.")
    
with col2_1:
    #st.image("static/images/Workflow.jpg")
    pass

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
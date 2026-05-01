import streamlit as st


st.set_page_config( page_title="AMP-PPI", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h1 style='text-align:center;'> AMP-PPI: A Machine Learning based tool for prediction of Protein-Protein Interaction between Antimicrobial Peptides and their Target Proteins</h1>", unsafe_allow_html=True)



st.markdown("---")
st.text("")

st.title("Tutorial")
with open("static/tutorial.pdf", "rb") as f:
    st.download_button("Download Tutorial File", f, "tutorial_amp-ppi.pdf")

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)

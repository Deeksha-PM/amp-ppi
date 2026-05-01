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
st.markdown(
    """
    <a href="static/tutorial.pdf" target="_blank" style="text-decoration: none;">
        <div style="
            background-color: #ff4b4b;
            color: white;
            padding: 20px;
            text-align: center;
            border-radius: 10px;
            font-size: 24px;
            font-weight: bold;
            margin: 20px 0;
            box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
        ">
            OPEN PDF TUTORIAL
        </div>
    </a>
    """,
    unsafe_allow_image_to_markdown=True,
    unsafe_allow_html=True
)

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)

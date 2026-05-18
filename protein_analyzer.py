import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from Bio.Seq import Seq
from Bio import SeqIO
import io

# st.markdown("""
#     <style>
#         /* Основной фон (тёмный режим) */
#         .stApp {
#             background: linear-gradient(135deg, #1a472a 0%, #0d2818 100%);
#         }
        
#         /* Карточки метрик */
#         [data-testid="stMetric"] {
#             background: rgba(255,255,255,0.1);
#             border-radius: 15px;
#             padding: 15px;
#             border: 1px solid rgba(255,255,255,0.2);
#             backdrop-filter: blur(10px);
#         }
        
#         /* Заголовки */
#         h1 {
#             color: #4ade80 !important;
#             text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
#             font-family: 'Courier New', monospace !important;
#         }
        
#         h2, h3 {
#             color: #86efac !important;
#         }
        
#         /* Текст в метриках */
#         [data-testid="stMetricValue"] {
#             color: #fbbf24 !important;
#             font-size: 2rem !important;
#         }
        
#         /* Кнопки */
#         .stButton > button {
#             background: linear-gradient(90deg, #10b981, #059669);
#             color: white;
#             border: none;
#             border-radius: 25px;
#             padding: 10px 25px;
#             font-weight: bold;
#             transition: all 0.3s;
#         }
        
#         .stButton > button:hover {
#             transform: scale(1.05);
#             box-shadow: 0 0 15px rgba(16,185,129,0.5);
#         }
        
#         /* Текстовое поле (аналог гелевой электрофореза) */
#         .stTextArea textarea {
#             background: #1e1e1e;
#             color: #10b981;
#             font-family: 'Monaco', 'Courier New', monospace;
#             border: 2px solid #10b981;
#             border-radius: 10px;
#         }
        
#         /* Прогресс-бары (для анализа) */
#         .stProgress > div > div {
#             background: linear-gradient(90deg, #10b981, #3b82f6);
#         }
#     </style>
# """, unsafe_allow_html=True)

st.title("🧬 Анализатор белковых последовательностей")
st.markdown("Введите последовательность белка в формате FASTA или просто вставьте текст:")

# Способ 1: текстовое поле
seq_input = st.text_area("Последовательность (без заголовка >)", height=150)

# Способ 2: загрузка файла
uploaded_file = st.file_uploader("Или загрузите FASTA файл", type=["fasta", "txt"])

if uploaded_file:
    content = uploaded_file.read().decode("utf-8")
    seq_input = "".join(content.splitlines()[1:])  # убираем header

def aa_frequency(seq):
    """Частота каждой аминокислоты (в % от длины)"""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    freq = {}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        freq[aa] = seq.count(aa) / len(seq) * 100
    return freq

def molecular_weight(seq):
    """Приблизительная молекулярная масса (Da)"""
    weights = {"A": 89, "R": 174, "N": 132, "D": 133, "C": 121,
               "Q": 146, "E": 147, "G": 75, "H": 155, "I": 131,
               "L": 131, "K": 146, "M": 149, "F": 165, "P": 115,
               "S": 105, "T": 119, "W": 204, "Y": 181, "V": 117}
    return sum(weights.get(aa, 0) for aa in seq.upper())

def charge_at_pH7(seq):
    """Заряд при pH 7: K+R - D-E"""
    seq = seq.upper()
    pos = seq.count("K") + seq.count("R")
    neg = seq.count("D") + seq.count("E")
    return pos - neg

if st.button("🔬 Анализировать"):
    if len(seq_input.strip()) == 0:
        st.warning("Введите последовательность")
    else:
        seq = seq_input.strip()
        
        # Частоты
        freq = aa_frequency(seq)
        df_freq = pd.DataFrame(freq.items(), columns=["Аминокислота", "Частота (%)"])
        
        # Метрики
        col1, col2, col3 = st.columns(3)
        col1.metric("Длина", f"{len(seq)} aa")
        col2.metric("Молекулярная масса", f"{molecular_weight(seq):.0f} Da")
        col3.metric("Заряд при pH 7", f"{charge_at_pH7(seq):+d}")
        
        # График частот
        fig = px.bar(df_freq, x="Аминокислота", y="Частота (%)", 
                     title="Аминокислотный состав")
        st.plotly_chart(fig)
        
        # Гидрофобность (гистограмма)
        hyd_dict = {"Hydrophobic": ["A","V","I","L","M","F","W","P"],
                    "Polar": ["S","T","Y","N","Q","C"],
                    "Charged": ["D","E","H","K","R"]}
        hyd_sums = {k: sum(freq[aa] for aa in v) for k,v in hyd_dict.items()}
        
        fig2, ax = plt.subplots()
        ax.pie(hyd_sums.values(), labels=hyd_sums.keys(), autopct="%1.1f%%")
        ax.set_title("Классы аминокислот")
        st.pyplot(fig2)
        
        # CSV экспорт
        st.download_button("📥 Скачать CSV", df_freq.to_csv(index=False), 
                           "aa_composition.csv", "text/csv")


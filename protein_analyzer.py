import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from Bio.Seq import Seq
from Bio import SeqIO
import io

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


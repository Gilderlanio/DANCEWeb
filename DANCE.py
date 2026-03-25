import streamlit as st
import pandas as pd
import networkx as nx
from st_cytoscape import cytoscape
from PIL import Image
from gprofiler import GProfiler

import matplotlib.cm as cm
import matplotlib.colors as mcolors



st.set_page_config(page_title = "Disease-Ancestry Networks", layout='wide')
logo = Image.open("imagens/logo.png")
st.sidebar.image(logo, use_container_width=True, width=100)

@st.cache_data
def load_data(path: str) -> pd.DataFrame:
    """Load and cache the GWAS dataset."""
    return pd.read_csv(path, sep="\t")

def get_pastel_red(value, vmin=0, vmax=1):
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.Reds  # escala de vermelho

    rgba = cmap(norm(value))

    # transformar em pastel (mistura com branco)
    pastel_factor = 0.5
    r, g, b, a = rgba
    r = r + (1 - r) * pastel_factor
    g = g + (1 - g) * pastel_factor
    b = b + (1 - b) * pastel_factor

    return mcolors.to_hex((r, g, b))

df = load_data("data/gwas_preprocessed_oct.tsv")
af_df = load_data("data/1000genomeFrequencies_alt.tsv")

st.sidebar.header("Filter")
select_risk = st.sidebar.checkbox("Filter for risk-alleles", value=True)
if select_risk:
    df = df[df["Odds"] > 1]

selected_phenotype = st.sidebar.multiselect("Select phenotype", sorted(df["Phenotype"].unique()))
if selected_phenotype:
    df = df[df["Phenotype"].isin(selected_phenotype)]

selected_author = st.sidebar.multiselect("Select author", sorted(df["Author"].unique()))
if selected_author:
    df = df[df["Author"].isin(selected_author)]

selected_gene = st.sidebar.multiselect("Select gene", sorted(df["Gene"].unique()))
if selected_gene:
    df = df[df["Gene"].isin(selected_gene)]

selected_region = st.sidebar.multiselect("Select region", sorted(df["Region"].unique()))
if selected_region:
    df = df[df["Region"].isin(selected_region)]

rsid_choice = st.sidebar.multiselect("Choose the rsID", options=[""] + sorted(df["rsID"].unique()))
if rsid_choice:
    df = df[df["rsID"].isin(rsid_choice)]

tab1, tab2, tab3, tab4 = st.tabs(["Data", "Network analysis", "Functional", "About"])

# GWAS data
with tab1:

    categorical_cols = ["Phenotype","Chr","Region","Gene","rsID"]
    categorical_summary = {}
    for col in categorical_cols:
        mode = df[col].mode()[0]
        freq = df[col].value_counts().max()
        unique = df[col].nunique()
        categorical_summary[col] = {"Most Common": mode, "Frequency": freq, "Unique": unique}

    # CSS customizado para headers
    st.markdown("""
    <style>
    
    button[data-baseweb="tab"] > div[data-testid="stMarkdownContainer"] p {
        font-size: 20px;  /* Change font size */
        font-weight: bold;  /* Optional */
        font-family: 'Georgia', serif;
    }
    
    .custom-header {
        font-family: 'Roboto', sans-serif;
        font-weight: 700;
        font-size: 28;
        color: #1b3a57;  /* azul moderno */
        margin-bottom: 0px;
    }
    
    .custom-subheader {
        font-family: 'Roboto', sans-serif;
        font-weight: 600;
        font-size: 24px;
        color: #1b3a57;  /* cinza escuro */
        margin-top: 20px;
        margin-bottom: 4px;
    }
    
    .custom-divider {
        border-top: 3px solid #1b3a57;
        margin: 10px 0 20px 0;
        width: 50px;
    }
    </style>
    """, unsafe_allow_html=True)

    # Aplicando estilos nos headers
    #st.markdown('<div class="custom-header">Disease-Ancestry Network (DANCE)</div>', unsafe_allow_html=True)
    #st.markdown('<div class="custom-divider"></div>', unsafe_allow_html=True)

    st.markdown('<div class="custom-subheader">Summary Statistics</div>', unsafe_allow_html=True)
    st.markdown('<div class="custom-divider"></div>', unsafe_allow_html=True)

    st.markdown("""
    <style>
    .nature-card {
        background-color: #fdfdfd;
        border-radius: 12px;
        padding: 10px;
        margin: 10px 5px;
        box-shadow: 0 3px 6px rgba(0,0,0,0.08);
        transition: transform 0.25s, box-shadow 0.25s;
        text-align: center;
        font-family: 'Georgia', serif;
    }
    .nature-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 10px 20px rgba(0,0,0,0.15);
    }
    .nature-card-title {
        font-weight: 600;
        font-size: 18px;
        color: #222222;
        margin-bottom: 4px;
    }
    .nature-card-value {
        font-size: 15px;
        font-weight: 700;
        color: #1b3a57; /* Tom azul escuro, elegante */
    }
    .nature-card-subvalue {
        font-size: 13px;
        font-weight: 500;
        color: #555555;
        margin-top: 0px;
    }
    </style>
    """, unsafe_allow_html=True)


    cards_per_row = 5
    cols_layout = st.columns(cards_per_row)

    st.markdown("---")

    for i, (col, stats) in enumerate(categorical_summary.items()):
        with cols_layout[i % cards_per_row]:
            st.markdown(f"""
            <div class="nature-card">
                <div class="nature-card-title">{col}</div>
                <div class="nature-card-value">{stats['Most Common']}</div>
                <div class="nature-card-subvalue">Frequency: {stats['Frequency']} | Unique: {stats['Unique']}</div>
            </div>
            """, unsafe_allow_html=True)

    st.dataframe(df, use_container_width=True, hide_index=True)

    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download", csv, "DANCE_filtered.csv", "text/csv")

#Network
with tab2:

    st.markdown(
        """
        <style>
        /* Seleciona todos os iframes dentro dos componentes */
        iframe {
            border: 3px solid #2E86C1;  /* cor e espessura da borda */
            border-radius: 15px;        /* cantos arredondados */
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);  /* sombra suave */
        }
        </style>
        """,
        unsafe_allow_html=True
    )
    if (selected_phenotype or selected_author or selected_gene or selected_region or rsid_choice):
        # usa agg_df para manter só o maior Odds por aresta

        populations = ["AFR", "AMR", "EAS", "EUR", "SAS"]
        pop = st.sidebar.selectbox("Population", populations, index=3)

        df["SNV"] = df["rsID"] + "-" + df["Allele"]
        af_df = af_df.rename(columns={"ID": "rsID"})
        af_df["SNV_ALT"] = af_df["rsID"] + "-" + af_df["ALT"]

        merged = df.merge(af_df, on="rsID", how="left", suffixes=("", "_af"))


        def get_correct_af(row):
            if pd.isna(row[pop]):
                return 0

            # risco = ALT
            if row["Allele"] == row["ALT"]:
                return row[pop]

            # risco = REF
            elif row["Allele"] == row["REF"]:
                return 1 - row[pop]

            # outro ALT não correspondente
            else:
                return 0

        merged["AF_corrected"] = merged.apply(get_correct_af, axis=1)

        merged["SNV"] = merged["rsID"] + "-" + merged["Allele"]
        af_map = merged.groupby("SNV")["AF_corrected"].max().to_dict()
        vmin = merged["AF_corrected"].min()
        vmax = merged["AF_corrected"].max()

        agg_df = df.groupby(["Phenotype", "SNV"], as_index=False)["Odds"].max()
        net = nx.from_pandas_edgelist(agg_df, source="Phenotype", target="SNV", edge_attr="Odds")

        # pega o maior Odds para normalização
        max_odds = agg_df["Odds"].max()

        # cor e espessura da aresta proporcional ao Odds
        for u, v, d in net.edges(data=True):
            d["color"] = "gray"
            d["width"] = 1 + 3 * (d["Odds"] / max_odds)

        # atributos dos nós
        for node in net.nodes():
            if node in df["Phenotype"].values:
                net.nodes[node]["color"] = "skyblue"
                net.nodes[node]["shape"] = "dot"
                net.nodes[node]["size"] = 50
            elif node in df["SNV"].values:
                af = af_map.get(node, 0)
                color = get_pastel_red(af, vmin, vmax)
                net.nodes[node]["color"] = color
                net.nodes[node]["size"] = 30
                net.nodes[node]["shape"] = "dot"

        elements = []
        # nós
        for node, data in net.nodes(data=True):
            elements.append({
                "data": {"id": str(node), "label": str(node)},
                "style": {
                    "background-color": data.get("color", "gray"),
                    "width": data.get("size", 20),
                    "height": data.get("size", 20),
                    "shape": data.get("shape", "ellipse"),
                },
            })
        # arestas
        for u, v, data in net.edges(data=True):
            elements.append({
                "data": {"source": str(u), "target": str(v)},
                "style": {
                    "line-color": data.get("color", "gray"),
                    "width": data.get("width", 2),
                },
            })
        # -------- Estilo base --------
        stylesheet = [
            {"selector": "node", "style": {"label": "data(label)", "background-color": "skyblue"}},
            {"selector": "edge", "style": {"curve-style": "bezier", "line-color": "gray"}}
        ]

        col1, col2 = st.columns(2)
        # Opções de layout

        layout_options = ["grid", "circle", "cose", "breadthfirst", "concentric"]
        layout = []
        with col1:
            selected_layout = st.radio("Layout:", layout_options, index=2, horizontal=True)

            cose_params = {
                "name": "cose",
                "fit": True,
                "padding": 10,
                "nodeRepulsion": 400000,
                "idealEdgeLength": 200,
                "edgeElasticity": 80,
                "gravity": 1,
                "numIter": 3000,
            }

            # botão reset
            if st.button("Reset"):
                if selected_layout == "cose":
                    layout = {**cose_params, "fit": True}
                else:
                    layout = {"name": selected_layout, "fit": True, "padding": 30}
            else:
                if selected_layout == "cose":
                    layout = {**cose_params, "fit": True}
                else:
                    layout = {"name": selected_layout, "fit": False}

        # Renderiza Cytoscape
        selected = cytoscape(elements, stylesheet, layout=layout, height="700px")

        st.markdown("**Node color = Risk Allele Frequency**")
        st.markdown("Light pink → Low frequency | Red → High frequency")



        if len(selected["nodes"]) > 0:
            selected_snps = [snp.split("-")[0] for snp in selected["nodes"] if snp.startswith("rs")]
            gp = GProfiler(return_dataframe=True)
            df_snp_sense = gp.snpense(selected_snps)
            st.dataframe(df_snp_sense)
    else:
        st.info(
            "The **DANCE** module renders the *SNP–Disease* network dynamically. "
            "The network is generated and displayed whenever a filter is selected from the sidebar."
        )

# Functional enrichment
with tab3:
    if len(df["Gene"].to_list()) <= 100:
        gp = GProfiler(return_dataframe=True)
        df_gene_enrichment = gp.profile(organism='hsapiens',
                   query=df["Gene"].to_list())

        target_cols =["source", "ID", "name", "p_value",
                      "significant", "description"]
        rename_map = {
            'source': 'Source',
            'native': 'ID',
            'name': 'Term',
            'term_size': 'Term Size',
            'query_size': 'Query Size',
            'p_value': 'P-Value',
            'significant': 'Significant?',
            'description': 'Description',
            'intersection_size': 'Overlap',
            'effective_domain_size': 'Domain Size',
            'precision': 'Precision',
            'recall': 'Recall',
            'query': 'Query Terms',
            'parents': 'Parents',
        }

        df_gene_enrichment = df_gene_enrichment.rename(columns=rename_map)
        st.dataframe(df_gene_enrichment, hide_index=True)
        df_snp_sense['variants'] = df_snp_sense['variants'].apply(
            lambda x: ", ".join(
                f"{k}: {v}"
                for k, v in x.items()
            ) if isinstance(x, dict) and x else ""
        )

        variant_labels = {
            "transcript_ablation": "Transcript ablation",
            "splice_acceptor_variant": "Splice acceptor",
            "splice_donor_variant": "Splice donor",
            "stop_gained": "Stop gained",
            "frameshift_variant": "Frameshift",
            "stop_lost": "Stop lost",
            "start_lost": "Start lost",
            "transcript_amplification": "Transcript amplification",
            "inframe_insertion": "Inframe insertion",
            "inframe_deletion": "Inframe deletion",
            "missense_variant": "Missense",
            "protein_altering_variant": "Protein altering",
            "splice_region_variant": "Splice region",
            "incomplete_terminal_codon_variant": "Incomplete terminal codon",
            "stop_retained_variant": "Stop retained",
            "synonymous_variant": "Synonymous",
            "coding_sequence_variant": "Coding sequence",
            "mature_miRNA_variant": "Mature miRNA",
            "5_prime_UTR_variant": "5' UTR",
            "3_prime_UTR_variant": "3' UTR",
            "non_coding_transcript_exon_variant": "Non-coding exon",
            "intron_variant": "Intron",
            "NMD_transcript_variant": "NMD transcript",
            "non_coding_transcript_variant": "Non-coding transcript",
            "upstream_gene_variant": "Upstream gene",
            "downstream_gene_variant": "Downstream gene",
            "TFBS_ablation": "TFBS ablation",
            "TFBS_amplification": "TFBS amplification",
            "TF_binding_site_variant": "TF binding site",
            "regulatory_region_ablation": "Regulatory region ablation",
            "regulatory_region_amplification": "Regulatory region amplification",
            "feature_elongation": "Feature elongation",
            "regulatory_region_variant": "Regulatory region",
            "feature_truncation": "Feature truncation",
            "intergenic_variant": "Intergenic",
            "splice_donor_5th_base_variant": "Splice donor +5 base",
            "splice_donor_region_variant": "Splice donor region",
            "splice_polypyrimidine_tract_variant": "Splice poly-pyrimidine tract"
        }

        # Converter contagens em flags (1 se > 0)
        flags = df_snp_sense['variants'].apply(lambda x: x.keys())

        flags_df = pd.json_normalize(flags).rename(columns=variant_labels)
        df_snp_sense = pd.concat([df_snp_sense, flags_df], axis=1)

        st.dataframe(df_snp_sense, hide_index=True)
    else:
        st.info(
            "The **DANCE** module renders the functional enrichment analysis dynamically."
            "This analysis is generated and displayed whenever a filter is selected from the sidebar."
        )

with tab4:
    st.header("Data sources")
    st.write("""
    The **DANCE** data is a subset of both the **1000 Genomes Project** data and the **NHGRI GWAS Catalog**.
    
    - **GWAS Catalog**: https://www.ebi.ac.uk/gwas/
    - **1000 Genomes Project**: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
    """)
    st.header("Citing DANCE")
    st.write("""
    If you use DANCE, please cite the following references:
    - **DANCE**: Araújo, GS, Lima, LH, Schneider, S, Leal, TP, da Silva, AP, Vaz de Melo, PO, Tarazona-Santos, E, Scliar, MO, Rodrigues, MR (2015). *Integrating, summarizing and visualizing GWAS-hits and human diversity with DANCE (Disease-ANCEstry Networks).* Bioinformatics. doi: [10.1093/bioinformatics/btv708](https://doi.org/10.1093/bioinformatics/btv708)
    - **1000 Genomes Project**: Siva, Nayanah. "1000 Genomes project." *Nature Biotechnology* 26.3 (2008): 256-256.
    - **NHGRI GWAS Catalog**: Welter, Danielle, et al. "The NHGRI GWAS Catalog, a curated resource of SNP-trait associations." *Nucleic Acids Research* 42.D1 (2014): D1001-D1006.
    """)

    st.header("Research Team")

    st.markdown("""
    
    - **Human Genetic Diversity Lab**  
      [LDGG](https://ldgh.com.br)
    
    - **Gilderlanio S. Araújo**  
      [ORCID](https://orcid.org/0000-0001-9199-9419) | [Lattes](http://lattes.cnpq.br/6152771446841901)
    - **Maíra R. Rodrigues**  
      [ORCID](https://orcid.org/0000-0003-3193-9558) | [Lattes](http://lattes.cnpq.br/1035025419267366)

    - **Eduardo Tarazona-Santos**  
      [ORCID](https://orcid.org/0000-0003-3508-3160) | [Lattes](http://lattes.cnpq.br/6203097295718656)
    """)

    st.header("Funding")

    capes = Image.open("imagens/capes.png")
    cnpq = Image.open("imagens/cnpq.png")
    fapespa = Image.open("imagens/fapespa.png")
    ufpa = Image.open("imagens/ufpa.png")
    ufmg = Image.open("imagens/ufmg.png")

    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.image(capes, use_container_width=False, width=70)

    with col2:
        st.image(cnpq, use_container_width=False, width=100)
    with col3:
        st.image(fapespa, use_container_width=False, width=100)
    with col4:
        st.image(ufmg, use_container_width=False, width=100)
    with col5:
        st.image(ufpa, use_container_width=False, width=100)

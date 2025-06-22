import streamlit as st
import pandas as pd
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
from rdkit import Chem
from rdkit.Chem import DataStructs, Draw
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit.Chem import Descriptors

Settings.CACHE = False

# Configuración de página
st.set_page_config(page_title="Drug Discovery Toolkit", layout="wide")

col1, col2 = st.columns([1, 6])
with col1:
    st.image("logo_app.png", width=650)
with col2:
    st.markdown("""
        <div style='background-color:#003366; padding:15px; border-radius:10px'>
            <h1 style='color:white; text-align:center;'>💊 Drug Discovery Toolkit</h1>
        </div>
    """, unsafe_allow_html=True)
    st.markdown("""
        <style>
        .stTabs [data-baseweb="tab"] {
            background-color: #e6f0ff;
            color: #003366;
            font-weight: bold;
        }
        .stTabs [aria-selected="true"] {
            background-color: #003366;
            color: white;
        }
        </style>
    """, unsafe_allow_html=True)

    st.markdown("Interfaz para explorar compuestos activos frente a dianas terapéuticas")

tabs = st.tabs([
    "IC50 (todos los compuestos)",
    "IC50 (compuestos preclínicos)",
    "Similitud con fármacos clínicos",
    "Visualización de estructura molecular"
])

# ---------------------- TAB 1: IC50 (todos los compuestos) ----------------------
with tabs[0]:
    st.header("Compuestos bioactivos frente una diana")

    target_query = st.text_input("🔎 Introduce el nombre o ChEMBL ID de una diana terapéutica (ej. EGFR o CHEMBL203):")
    ic50_threshold = st.slider("🎚️ Filtrar compuestos con IC50 menor a (nM):", 0, 1000, 100, step=1)
    max_to_show = st.slider("🔢 Número de compuestos a mostrar (ordenados por menor IC50):", 10, 1000, 100, step=10)

    if target_query:
        with st.spinner("🔍 Buscando diana y compuestos..."):
            if target_query.upper().startswith("CHEMBL"):
                targets = [new_client.target.get(target_query.upper())]
            else:
                targets = list(new_client.target.filter(target_synonym__icontains=target_query)) + \
                          list(new_client.target.filter(pref_name__icontains=target_query))

        if targets:
            target_names = [f"{t['pref_name']} ({t['target_chembl_id']})" for t in targets]
            selected = st.selectbox("🎯 Selecciona la diana:", target_names)
            target_id = selected.split('(')[-1][:-1]

            with st.spinner("📡 Descargando bioactividades IC50..."):
                activities = new_client.activity.filter(target_chembl_id=target_id).only(
                    ["molecule_chembl_id", "standard_value", "standard_type", "standard_units", "canonical_smiles"]
                )

                data = []
                seen_ids = set()

                for a in activities:
                    mol_id = a.get("molecule_chembl_id")
                    tipo = a.get("standard_type", "").lower()
                    units = a.get("standard_units")
                    value_raw = a.get("standard_value")
                    smiles = a.get("canonical_smiles")

                    if not mol_id or mol_id in seen_ids or tipo != "ic50" or not value_raw or not smiles:
                        continue

                    try:
                        value = float(value_raw)
                        if units == "nM":
                            value_nM = value
                        elif units == "µM":
                            value_nM = value * 1000
                        elif units == "pM":
                            value_nM = value / 1000
                        else:
                            continue
                        if not Chem.MolFromSmiles(smiles):
                            continue
                    except:
                        continue

                    data.append({"ChEMBL_ID": mol_id, "IC50 (nM)": value_nM, "SMILES": smiles})
                    seen_ids.add(mol_id)

                df = pd.DataFrame(data)
                df = df[df["IC50 (nM)"] <= ic50_threshold].sort_values("IC50 (nM)").reset_index(drop=True)
                df_top = df.head(max_to_show)

                if not df_top.empty:
                    st.success(f"✅ {len(df_top)} compuestos encontrados")
                    st.dataframe(df_top, use_container_width=True)

                    ref = df_top.iloc[0]
                    st.markdown("### 🔬 Compuesto con menor IC50")
                    st.code(f"ID: {ref['ChEMBL_ID']}\nIC50: {ref['IC50 (nM)']} nM\nSMILES: {ref['SMILES']}")
                    mol_img = Chem.MolFromSmiles(ref['SMILES'])
                    if mol_img:
                        st.image(Draw.MolToImage(mol_img, size=(300, 300)), caption="Estructura molecular")

                    csv = df_top.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Descargar CSV", data=csv, file_name="resultados_ic50.csv", mime="text/csv")

                    # ---------------------- ANÁLISIS EXPLORATORIO  ----------------------
                    st.markdown("## 📊 Estadísticas del IC50")

                    import io

                    fig_ic50 = plt.figure(figsize=(5, 3))
                    sns.histplot(df_top['IC50 (nM)'], kde=True, bins=30)
                    plt.xlabel("IC50 (nM)")
                    plt.ylabel("Frecuencia")
                    plt.title("Distribución de IC50")
                    st.pyplot(fig_ic50, use_container_width=False)

                    # Botón para descargar
                    buf = io.BytesIO()
                    fig_ic50.savefig(buf, format="png", bbox_inches='tight')
                    st.download_button(
                        label="📥 Descargar gráfico IC50",
                        data=buf.getvalue(),
                        file_name="ic50_distribution.png",
                        mime="image/png",
                        key = "descargar_grafico_ic50_tab0"
                    )



                else:
                    st.warning("⚠️ No se encontraron compuestos con IC50 dentro del umbral.")
        else:
            st.error("❌ No se encontró ninguna diana con ese nombre o ID.")

# ---------------------- TAB 2: IC50 (solo no clínicos) ----------------------
with tabs[1]:
    st.header("Compuestos preclínicos bioactivos frente una diana")

    target_query_nc = st.text_input("🔎 Introduce el nombre o ChEMBL ID de una diana terapéutica (ej. EGFR o CHEMBL203):",
                                    key="preclin_target")
    ic50_threshold_nc = st.slider("🎯 Filtrar compuestos con IC50 menor a (nM):", 1, 1000, 100, step=10,
                                  key="preclin_ic50")
    max_to_show_nc = st.slider("🔢 Número de compuestos a mostrar (ordenados por menor IC50):", 10, 1000, 100,
                               step=10, key="preclin_max")

    if target_query_nc:
        with st.spinner("🔍 Buscando diana y compuestos..."):
            if target_query_nc.upper().startswith("CHEMBL"):
                targets = [new_client.target.get(target_query_nc.upper())]
            else:
                targets = list(new_client.target.filter(target_synonym__icontains=target_query_nc)) + \
                          list(new_client.target.filter(pref_name__icontains=target_query_nc))

        if not targets:
            st.error("❌ No se encontró ninguna diana con ese nombre o ID.")

        if targets:  # Asegurarse de que targets no esté vacío antes de continuar
            target_names = [f"{t['pref_name']} ({t['target_chembl_id']})" for t in targets]
            selected = st.selectbox("🎯 Selecciona la diana encontrada:", target_names, key="preclin_select")
            target_id = selected.split('(')[-1][:-1]

            with st.spinner("📡 Descargando actividades IC50 desde ChEMBL..."):
                activities = new_client.activity.filter(target_chembl_id=target_id).only(
                    ["molecule_chembl_id", "standard_value", "standard_type", "standard_units", "canonical_smiles"]
                )

                data = []
                seen_ids = set()

                for a in activities:
                    mol_id = a.get("molecule_chembl_id")
                    tipo = a.get("standard_type", "").lower()
                    units = a.get("standard_units")
                    value_raw = a.get("standard_value")
                    smiles = a.get("canonical_smiles")

                    if not mol_id or mol_id in seen_ids or tipo != "ic50" or not value_raw or not smiles:
                        continue

                    try:
                        value = float(value_raw)
                        if units == "nM":
                            ic50 = value
                        elif units == "µM":
                            ic50 = value * 1000
                        elif units == "pM":
                            ic50 = value / 1000
                        else:
                            continue

                        if not Chem.MolFromSmiles(smiles):
                            continue

                        data.append({"ChEMBL_ID": mol_id, "IC50 (nM)": ic50, "SMILES": smiles})
                        seen_ids.add(mol_id)

                    except Exception as e:
                        continue

            if not data:
                st.warning("⚠️ No se encontraron actividades IC50 válidas para esta diana.")
                # No st.stop() aquí
            else:
                df = pd.DataFrame(data)
                df = df[df["IC50 (nM)"] <= ic50_threshold_nc]

                # Cachear la función para evitar llamadas repetidas a la API
                @st.cache_data(ttl=3600)  # Cachear por 1 hora
                def get_non_clinical_ids(chembl_ids):
                    non_clin_ids = set()
                    batch_size = 500
                    for i in range(0, len(chembl_ids), batch_size):
                        batch = chembl_ids[i:i + batch_size]
                        try:
                            # Se asegura que el objeto 'molecule' se inicialice correctamente en memoria (caché)
                            # Se instancia 'new_client.molecule' dentro de la función para garantizar su correcta serialización en caché
                            mol_api = new_client.molecule
                            query = mol_api.filter(molecule_chembl_id__in=batch).only("molecule_chembl_id", "max_phase")
                            for mol in query:
                                phase = mol.get("max_phase")
                                # Considera un compuesto como preclínico si max_phase es 0 o None (no ha pasado a ensayos clínicos)
                                if phase in [None, 0, '0', 'Preclinical', 'preclinical']:
                                    non_clin_ids.add(mol["molecule_chembl_id"])
                        except Exception as e:
                            st.warning(
                                f"Error al consultar la API de ChEMBL para fase clínica: {e}. Algunos compuestos podrían no ser filtrados correctamente.")
                            continue  # Permite continuar con el siguiente lote si ocurre un error
                    return non_clin_ids


                if not df.empty:
                    ids_to_check = df["ChEMBL_ID"].dropna().unique().tolist()
                    with st.spinner("⏳ Verificando fase clínica de los compuestos..."):
                        preclin_ids = get_non_clinical_ids(ids_to_check)

                    df = df[df["ChEMBL_ID"].isin(preclin_ids)].copy()

                if df.empty:
                    st.info("⚠️ No se encontraron compuestos preclínicos con los criterios seleccionados.")
                else:
                    df = df.sort_values("IC50 (nM)").reset_index(drop=True)
                    df_top = df.head(max_to_show_nc)

                    st.success(f"✅ {len(df_top)} compuestos preclínicos mostrados")
                    st.dataframe(df_top, use_container_width=True)

                    if not df_top.empty:
                        ref = df_top.iloc[0]
                        st.markdown("### 🔬 Compuesto preclínico con menor IC50")
                        st.code(f"ID: {ref['ChEMBL_ID']}\nIC50: {ref['IC50 (nM)']} nM\nSMILES: {ref['SMILES']}")
                        mol_img = Chem.MolFromSmiles(ref['SMILES'])
                        if mol_img:
                            st.image(Draw.MolToImage(mol_img, size=(300, 300)), caption="Estructura molecular")
                        else:
                            st.warning(
                                "No se pudo generar la imagen de la estructura molecular para el compuesto más potente (SMILES inválido).")

                    csv = df_top.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Descargar CSV", csv, file_name="compuestos_preclinicos.csv", mime="text/csv")

                    # ---------------------- ANÁLISIS EXPLORATORIO  ----------------------
                    st.markdown("## 📊 Estadísticas del IC50")

                    import io

                    fig_ic50 = plt.figure(figsize=(5, 3))
                    sns.histplot(df_top['IC50 (nM)'], kde=True, bins=30)
                    plt.xlabel("IC50 (nM)")
                    plt.ylabel("Frecuencia")
                    plt.title("Distribución de IC50")
                    st.pyplot(fig_ic50, use_container_width=False)

                    # Botón para descargar
                    buf = io.BytesIO()
                    fig_ic50.savefig(buf, format="png", bbox_inches='tight')
                    st.download_button(
                        label="📥 Descargar gráfico IC50",
                        data=buf.getvalue(),
                        file_name="ic50_distribution.png",
                        mime="image/png",
                        key = "descargar_grafico_ic50_tab1"
                    )




    else:  # Si targets está vacío (no se encontró la diana)
            st.error("❌ No se encontró ninguna diana con ese nombre o ID.")

# ---------------------- TAB 3: Similitud con fármacos clínicos ----------------------
with tabs[2]:
    st.header("Fármacos clínicos similares para el reposicionamiento")

    modo = st.radio("Modo de entrada de compuesto de referencia:", ["SMILES manual", "Subir CSV"])

    @st.cache_data(show_spinner=True, ttl=86400)
    def get_all_clinical_molecules(max_mols_to_fetch):
        clinical = []
        query_iterator = new_client.molecule.filter(
            max_phase__in=[1, 2, 3, 4],
            molecule_structures__isnull=False
        ).only(["molecule_chembl_id", "pref_name", "max_phase", "molecule_structures"])

        fetched_count = 0
        for mol in query_iterator:
            if fetched_count >= max_mols_to_fetch:
                break

            smiles = mol.get("molecule_structures", {}).get("canonical_smiles")
            chembl_id = mol["molecule_chembl_id"]

            if smiles:
                mol_rdkit = Chem.MolFromSmiles(smiles)
                if mol_rdkit:
                    indication = "No especificada"
                    try:
                        indications = new_client.drug_indication.filter(molecule_chembl_id=chembl_id)
                        if indications:
                            indication = "; ".join(set([
                                ind.get("mesh_heading") or ind.get("efo_term") or "N/A"
                                for ind in indications if ind.get("mesh_heading") or ind.get("efo_term")
                            ])) or "No disponible"
                    except Exception:
                        indication = "Error al obtener indicación"

                    clinical.append({
                        "ChEMBL_ID": chembl_id,
                        "Nombre": mol.get("pref_name") or "(Sin nombre)",
                        "Fase_Clinica": mol.get("max_phase", "Desconocido"),
                        "Indicación": indication,
                        "SMILES": smiles,
                        "Mol": mol_rdkit
                    })
                    fetched_count += 1
        return clinical

    max_mols = st.slider("🔢 Número máximo de fármacos clínicos a comparar:", 100, 10000, 1000, step=100)

    with st.spinner("⏳ Cargando fármacos clínicos..."):
        clinical_mols = get_all_clinical_molecules(max_mols)
        for mol in clinical_mols:
            if mol["Mol"]:
                mol["Fingerprint"] = GetMorganFingerprintAsBitVect(mol["Mol"], 2, nBits=2048)
            else:
                mol["Fingerprint"] = None

    if modo == "SMILES manual":
        input_smiles = st.text_input("Introduce el SMILES del compuesto de interés:")
        if input_smiles:
            mol = Chem.MolFromSmiles(input_smiles)
            if not mol:
                st.error("❌ SMILES inválido.")
            else:
                fp_query = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                resultados = []
                with st.spinner("⏳ Calculando similitud..."):
                    for target in clinical_mols:
                        if target.get("Fingerprint"):
                            sim = DataStructs.TanimotoSimilarity(fp_query, target["Fingerprint"])
                            resultados.append({
                                "SMILES_Entrada": input_smiles,
                                "ChEMBL_ID_Clinico": target["ChEMBL_ID"],
                                "Nombre_Clinico": target["Nombre"],
                                "Fase_Clinica": target["Fase_Clinica"],
                                "Indicación": target.get("Indicación"),
                                "SMILES_Clinico": target["SMILES"],
                                "Similaridad_Tanimoto": sim
                            })

                if resultados:
                    df = pd.DataFrame(resultados).sort_values("Similaridad_Tanimoto", ascending=False).reset_index(drop=True)
                    sim_threshold = st.slider("Filtrar resultados con similitud mínima:", 0.0, 1.0, 0.3, step=0.05)
                    df = df[df["Similaridad_Tanimoto"] >= sim_threshold]
                    st.success(f"🔍 {len(df)} resultados tras filtrar")
                    st.dataframe(df.head(100), use_container_width=True)

                    csv_out = df.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Descargar CSV", csv_out, file_name="similitud_resultados.csv", mime="text/csv")
                else:
                    st.info("No hay compuestos clínicos similares que cumplan el umbral.")
    else:
        file = st.file_uploader("📁 Sube un CSV con columna 'SMILES'", type=["csv"])
        if file:
            df_input = pd.read_csv(file)
            if "SMILES" not in df_input.columns:
                st.error("❌ El archivo debe tener una columna llamada 'SMILES'.")
            else:
                resultados = []
                with st.spinner("⏳ Comparando compuestos..."):
                    for i, row in df_input.iterrows():
                        query_smiles = row["SMILES"]
                        mol = Chem.MolFromSmiles(query_smiles)
                        id_query = row.get("ChEMBL_ID", f"Mol_{i+1}")
                        if mol:
                            fp_query = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                            for target in clinical_mols:
                                if target.get("Fingerprint"):
                                    sim = DataStructs.TanimotoSimilarity(fp_query, target["Fingerprint"])
                                    resultados.append({
                                        "ChEMBL_ID_Entrada": id_query,
                                        "SMILES_Entrada": query_smiles,
                                        "ChEMBL_ID_Clinico": target["ChEMBL_ID"],
                                        "Nombre_Clinico": target["Nombre"],
                                        "Fase_Clinica": target["Fase_Clinica"],
                                        "Indicación": target.get("Indicación"),
                                        "SMILES_Clinico": target["SMILES"],
                                        "Similaridad_Tanimoto": sim
                                    })

                if resultados:
                    df = pd.DataFrame(resultados).sort_values("Similaridad_Tanimoto", ascending=False).reset_index(drop=True)
                    st.success(f"🔍 Comparaciones completadas: {len(df)} pares")
                    st.dataframe(df.head(100), use_container_width=True)
                    csv_out = df.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Descargar CSV", csv_out, file_name="similitud_resultados.csv", mime="text/csv")

# ---------------------- TAB 4: Visualización de estructuras ----------------------

import io

with tabs[3]:
    st.header("Comparar estructuras moleculares a partir de SMILES")

    col_inputs = st.columns(2)
    with col_inputs[0]:
        smiles1 = st.text_input("🔹 Introduce el primer SMILES:", key="smiles1")
    with col_inputs[1]:
        smiles2 = st.text_input("🔸 Introduce el segundo SMILES:", key="smiles2")

    if smiles1 or smiles2:
        col_mols = st.columns(2)

        if smiles1:
            mol1 = Chem.MolFromSmiles(smiles1)
            with col_mols[0]:
                if mol1:
                    st.subheader("🧪 Estructura 1")
                    img1 = Draw.MolToImage(mol1, size=(300, 300))
                    buf1 = io.BytesIO()
                    img1.save(buf1, format="PNG")
                    st.image(img1, caption="Molécula 1")
                    st.download_button("📥 Descargar estructura 1", data=buf1.getvalue(),
                                       file_name="mol1.png", mime="image/png", key="descargar1")
                else:
                    st.error("❌ El primer SMILES no es válido.")

        if smiles2:
            mol2 = Chem.MolFromSmiles(smiles2)
            with col_mols[1]:
                if mol2:
                    st.subheader("🧪 Estructura 2")
                    img2 = Draw.MolToImage(mol2, size=(300, 300))
                    buf2 = io.BytesIO()
                    img2.save(buf2, format="PNG")
                    st.image(img2, caption="Molécula 2")
                    st.download_button("📥 Descargar estructura 2", data=buf2.getvalue(),
                                       file_name="mol2.png", mime="image/png", key="descargar2")
                else:
                    st.error("❌ El segundo SMILES no es válido.")

#pie de pagina
st.markdown("""
    <hr style='border-top: 1px solid #bbb; margin-top: 30px;'>
    <p style='text-align: center; color: grey; font-size: small;'>
        Aplicación desarrollada como parte del Trabajo de Fin de Grado – Universidad Pablo de Olavide, 2025
    </p>
""", unsafe_allow_html=True)

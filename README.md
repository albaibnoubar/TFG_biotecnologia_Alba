# TFG_biotecnologia_Alba

**Aplicación para el análisis computacional de fármacos contra proteínas diana**

Este repositorio contiene la aplicación desarrollada como parte del Trabajo de Fin de Grado en el grado de Biotecnología. La herramienta permite consultar datos de bioactividad (IC50) de compuestos frente a proteínas diana, filtrar por fase clínica y visualizar resultados relevantes desde la base de datos ChEMBL.

## Funcionalidades

- Búsqueda de proteínas diana por nombre o ID ChEMBL
- Descarga de datos IC50 desde ChEMBL mediante su API
- Filtrado por fase clínica 
- Visualización y exportación de compuestos activos en formato CSV
- Similitud con fármacos clínicos basada en estructuras moleculares

## Tecnologías usadas

- Python 3.12
- Streamlit
- RDKit
- pandas
- ChEMBL Web Resource Client

## Estructura del repositorio
```
CODIGO/
├── DDT.py # Código principal de la app
├── logo_app.png # Logotipo de la aplicación
├── requirements.txt # Librerías necesarias para ejecutar la app
└── README.md # Descripción del proyecto
```

## Cómo ejecutar la aplicación

1. Clona el repositorio:
   ```bash
   git clone https://github.com/albainoubar/TFG_biotecnologia_Alba.git
   ```
   
2. Instala los requisitos
   ```bash
   pip install -r requirements.txt
   ```

3. Ejecuta la app con Streamlit
   ```bash
   streamlit run DDT.py
   ```


# 🧬 Gene Expression & Survival Analysis in Prostate Cancer

An interactive **R Shiny web application** for exploring bulk RNA-Seq gene expression and evaluating the **prognostic significance of genes in prostate adenocarcinoma (PRAD)** using survival analysis.

The app is built using **TCGA-PRAD data** and enables **no-code biomarker exploration** for researchers, clinicians, and bioinformaticians.

---

## 🌐 Live Application

👉 **Access the deployed app here:**  
🔗 https://sgrbnf.shinyapps.io/PRAD_SURV/

No installation required — runs directly in your browser.

---

## 🚀 Key Features

- 📊 **Gene expression visualization**
  - Tumor vs normal comparison using boxplots
- ⏱ **Survival analysis**
  - Kaplan–Meier curves
  - Median-based stratification (UP vs DOWN expression)
  - Log-rank test p-values
- 📂 **Transparent & reproducible**
  - All datasets used in the app are provided in this repository
- 📥 **Downloadable results**
  - Expression and survival plots exported as PDF
- 🎨 **Modern UI**
  - Built with `bslib` (Bootstrap 5)

---

## 📁 Repository Contents

```
├── app.R
├── README.md
├── clinical_data_surv.txt
├── TopGenes_HUGO.txt
├── v_normsd.txt
├── www/
│   └── wrk2.jpg
└── LICENSE
```
---

## 🛠 Run the App Locally (Optional)

To install nad run the nextflow pipeline, follow these steps:

1. Clone this repository:

```bash
git clone https://github.com/<your-username>/PRAD_SURV.git
cd PRAD_SURV
```

2. Install required R packages:

```bash
install.packages(c(
  "shiny",
  "dplyr",
  "survival",
  "survminer",
  "ggplot2",
  "bslib"
))

```
3. Run the app:
```bash
shinyApp(ui = ui, server = server)
```
## Workflow:

<img src="https://github.com/sgr308/PRAD_SURV/blob/main/www/wrk2.jpg?raw=true"/>

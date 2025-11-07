This repository provides the complete R code used to generate the main and supplementary figures for the study:
“Longitudinal Genetic Architecture of Bone Mineral Density Trajectories and Fracture Risk in Postmenopausal Women.”

Overview

The project integrates multi-cohort genomic and phenotypic analyses to identify genetic determinants of longitudinal bone mineral density (BMD) trajectories and their impact on fracture risk in postmenopausal women. The repository includes visualization scripts for Manhattan plots, Mendelian randomization (MR) forest plots, and functional single-cell RNA-seq validation.

Included Figures

Figure 1: Manhattan plot for Hip BMD (Black population, WHI TrajGWAS)

Figure 2: Manhattan plot for Spine BMD (Black population, WHI TrajGWAS)

Figure 3: SPTBN1 single-cell expression across human tissues (GSE134355)

Figure 4: Functional GO enrichment and developmental pseudotime analysis of SPTBN1

Figure 5: One-sample Mendelian randomization of genetically predicted BMD and fracture risk (FHS)

Data Sources

Women’s Health Initiative (WHI): phs000200.v13.p3

Framingham Heart Study (FHS): phs000007.v32.p13

Single-cell RNA-seq dataset (Adult Human Cell Atlas): GSE159929 and GSE134355

All raw data are available via their respective repositories (dbGaP, GEO). This Zenodo deposit provides only scripts and derived figure outputs, not individual-level data.

Software and Requirements

R (≥4.2.0)

Packages: ggplot2, dplyr, patchwork, cowplot, ComplexHeatmap, Seurat, ggrepel, and dependencies.
Each figure script is self-contained and saves high-resolution PNG and PDF outputs to the output_figures/ directory.

Ethics and Data Access

All analyses were performed under IRB approval (Ohio State University, STUDY20250826) with WHI and FHS Data Use Agreements. Public GEO data (GSE134355) were used for transcriptomic validation.

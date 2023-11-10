# SomaLogic_Benchmark

This repository contains all custom scripts created during the SomaLogic Benchmark project as part of the University of Oregon Bioinformatics and Genomics Masters Program (UO BGMP) internship at Oregon Health & Science University Cancer Early Detection Advanced Research Center (OHSU CEDAR). Data is not included.

**Project Goal:** To benchmark the performance of SomaLogic SomaScan (7596 aptamer panel) against Seer Proteograph MS for proteomic analysis of murine and human data.

**Background & Results:** See [2023_Fall_Soriano](2023_Fall_Soriano.docx) for project details, methods, and results.

**Repository Guide:**

| File Name | Description |
| :---:   | :---: |
| [Somalogic_Benchmark_v4.Rmd](Somalogic_Benchmark_v4.Rmd) | R markdown for main project analysis. |
| [NomPivot_MouseHuman.R](NomPivot_MouseHuman.R) | R script with custom functions for mouse-human nomenclature pivots based on submitted gene names. Not currently used in main script, but can be used if analysis is targeted at gene-based results instead of protein-based results. |
| [NomPivot_MouseHuman_UniProtInput.R](NomPivot_MouseHuman_UniProtInput.R) | R script with custom functions for mouse-human nomenclature pivots based on submitted UniProt IDs. Currently used in main script for protein-based analysis. |
| [Seer_HumanCC_WiT.R](Seer_HumanCC_WiT.R) | R script used to process human Seer FC data from main OHSU CEDAR prostate cancer proteomic study. |
| [LabNotebook_SomaBenchmark_SS.txt](LabNotebook_SomaBenchmark_SS.txt) | Lab notebook for SomaLogic Benchmark Project. |
| [2023_Fall_Soriano.docx](2023_Fall_Soriano.docx) | Fall term paper submitted for UO BGMP - contains project background and discussion on relevant results. |

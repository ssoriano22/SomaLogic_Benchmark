# SomaLogic_Benchmark - 2024 Updates (Milestone 2 Work)

This repository contains all updated scripts, an updated lab notebook, and final Milestone 2 SomaLogic Benchmark presentation. Data is not included.

**Project Goal:** To benchmark the performance of SomaLogic SomaScan (7596 aptamer panel) against Seer Proteograph MS for proteomic analysis of murine and human data.

**Background & Results:** See [2023_Fall_Soriano](2023_Fall_Soriano.docx) for project details, methods, and results from first Milestone/UO BGMP internship work.

**Repository Guide:**

| File Name | Description |
| :---:   | :---: |
| [Somalogic_Benchmark_v6.Rmd](Somalogic_Benchmark_v6.Rmd) | R markdown for main project analysis. Uses [NomPivot_MouseHuman_UniProtInput.R](../NomPivot_MouseHuman_UniProtInput.R) from earlier work on this project. |
| [Somalogic_Benchmark_v6.html](Somalogic_Benchmark_v6.html) | Knitted R markdown for main project analysis. Simple captions are provided for each figure/set of related figures. |
| [Seer_HumanCC_WiT.R](Seer_HumanCC_WiT.R) | R script used to process human Seer data from main OHSU CEDAR prostate cancer proteomic study. Revised to match murine data processing pipeline and include MaxQuant initial processing steps. |
| [SB_RawSeerDataProcessing.Rmd](SB_RawSeerDataProcessing.Rmd) | R markdown used to process murine Seer data from main OHSU CEDAR murine PDAC project. Revised to match human data processing pipeline. |
| [LabNotebook_SomaBenchmarkcon_SS.txt](LabNotebook_SomaBenchmarkcon_SS.txt) | Lab notebook for SomaLogic Benchmark Project. |
| [SomaLogic_Milestone2_20240209_v4.pptx](SomaLogic_Milestone2_20240209_v4.pptx) | Final version of Milestone 2 presentation - includes some figures created by Dr. Matthew Chang using this correlation data, but the plotting code is not located here.|

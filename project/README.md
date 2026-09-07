# Guideline-Aware LLM Extraction for Gastric Cancer Risk Stratification

## Overview
This repository contains the code, synthetic v5 dataset, fine-tuned model, evaluation and output results for the thesis *"Guideline-Aware LLM Extraction for Gastric Cancer Risk Stratification: Development and Synthetic Validation"* (Zhang et al., 2026).

This project presents an AI pipeline for automated gastric cancer risk stratification following the **MAPS III (2025)** clinical guideline. The system:

- Extracts **39 structured clinical features** from free‑text endoscopy and pathology reports using a **fine‑tuned Qwen2.5‑1.5B** large language model (with QLoRA).
- Maps extracted features to guideline‑based management recommendations via a **deterministic decision tree**.
- Provides **traceable outputs** with explicit MAPS III reference numbers (e.g., REC.31) and warnings.

All training and evaluation were performed on a **synthetically generated corpus (V5)** that adheres to MAPS III reporting standards. The pipeline is **open‑source** and designed for local, privacy‑preserving deployment.

---
## Repository Structure

```
.
├── README.md                         
├── LICENSE                           
├── requirements.txt                  
├── .gitignore                        
│
├── data/                             # V5 dataset
│   └── raw_reports&ground_truth.zip                #  200 patients（.txt）and json features with values
│ 
├── patient_report_generator_code/ 
│   └── apiprompt2.ipynb
│
├── schema/                           
│   └── extraction_schema_2025_sc_39new.json        
│
├── alias_mapping/                          # Output normalization mapping
│   └── qwen_output_alias_map4_1.json
│
├── models/                           
│   └── lora_adapter.zip
│
├── decision_tree_mapping/                           
│   └── mapping.py
│
├──notebooks/                        # Core code
│   ├── 01_finetune_qwen_extraction.ipynb       # Model fine-tuning, clinical field extraction, and accuracy evaluation
│   └── 02_integrate.ipynb   # Integrate the pipeline and generate the final recommendation report
└── example_output/
    ├── extraction/
    └── recommendation/

```
## Quick Start

All code for this project is provided in **Jupyter Notebook** (`.ipynb`) format. Please follow the steps below:

### Core Execution Workflow
Perform the following steps according to different requirements:

1. **`01_finetune_qwen_extraction.ipynb`**
    - Fine-tune the model from scratch
    - Output structured JSON extraction results 
    - Output extraction accuracy evaluation

2. **`02_integrate.ipynb`**
   - Directly use the fine-tuned model for patient classification

> 💡 **Tip**: If only need to view the final results, directly open `example_output` and skip the training and classification steps

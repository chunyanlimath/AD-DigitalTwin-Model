# Data-driven spatiotemporal modeling reveals personalized trajectories of cortical atrophy in Alzheimer's disease

This repository contains the source code, trained models, and supplementary scripts for the manuscript:

> **Data-driven spatiotemporal modeling reveals personalized trajectories of cortical atrophy in Alzheimer's disease**

---

## 🔍 Overview
This project develops a **spatiotemporal network-based PDE model** to describe and predict the progression of **Alzheimer’s Disease (AD)** biomarkers — **Aβ**, **tau**, **neurodegeneration**, and **cognitive decline** — using longitudinal **neuroimaging data** and **brain functional connectivity network**.

The model formulates the biomarker cascade as a system of coupled PDEs defined on the **68-region Desikan–Killiany brain atlas functional connectivity network**. It captures both spatial propagation and temporal evolution of biomarker dynamics, enabling personalized forecasting and regional sensitivity analysis.

---

## 🧠 Model Highlights
- **Network-based PDE formulation:**  
  \[
  \frac{dy}{dt} + D L^G y = f(y; p)
  \]
  where \(L^G\) is the graph Laplacian derived from the subject-specific functional connectivity matrix.
  
- **Hierarchical Training with Homotopy Regularization:**  
  Efficiently calibrates model parameters using neuroimaging data (Aβ-PET, tau-PET, cortical thickness, MMSE).

- **Two-level Sobol Sensitivity Analysis:**  
  Identifies **critical brain regions and parameter interactions** driving disease progression.

- **Personalized Digital Twin:**  
  Enables individualized biomarker forecasting and identification of **sensitive targets for potential therapy**.

---

## 📊 Data
The model is trained and validated on data from the **Alzheimer’s Disease Neuroimaging Initiative (ADNI)**.  
To access ADNI data, please register at: [https://adni.loni.usc.edu](https://adni.loni.usc.edu)

---

## 🧩 Repository Structure

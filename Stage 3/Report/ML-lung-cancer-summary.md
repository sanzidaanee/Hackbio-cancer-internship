# 🫁 Machine Learning for Lung Cancer Biomarker Identification

**Project:** Identify Potential Biomarkers for Cancer Detection using ML Models  
**Dataset:** TCGA Lung Adenocarcinoma (LUAD)  
**Language:** Python  
**Notebook:** `ML_Lung_Cancer.ipynb`

---

## 📌 Overview

This module applies supervised machine learning to classify lung tissue samples as **Primary Tumor** or **Solid Tissue Normal** using gene expression data from The Cancer Genome Atlas (TCGA) LUAD dataset. By combining **Recursive Feature Elimination (RFE)** for feature selection with a **Random Forest classifier**, the pipeline identifies the most discriminating gene expression features — potential biomarkers for lung adenocarcinoma detection.

---

## 🎯 Objectives

- Predict sample type (Primary Tumor vs. Solid Tissue Normal) from gene expression profiles
- Reduce feature dimensionality using RFE to identify the most informative genes
- Train and evaluate a Random Forest classifier on the selected features
- Identify candidate biomarker genes relevant to lung cancer diagnosis

---

## 📂 Dataset

| Property | Detail |
|----------|--------|
| Source | The Cancer Genome Atlas (TCGA) |
| Cancer Type | Lung Adenocarcinoma (LUAD) |
| Total Samples | 521 (subset of 40 used: 20 Primary Tumor + 20 Solid Tissue Normal) |
| Features | Normalized gene expression values (RNA-Seq counts) |
| Target | Sample type: `Primary Tumor` / `Solid Tissue Normal` |

---

## 🔬 Workflow

```
<img width="131" height="150" alt="ml_lung_cancer_workflow" src="https://github.com/user-attachments/assets/11fb28a2-2c5d-436a-b0b8-8c755fe00024" />




<svg width="100%" viewBox="0 0 680 780" role="img" style="" xmlns="http://www.w3.org/2000/svg">
  <title style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">ML Lung Cancer Classification Pipeline</title>
  <desc style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">Flowchart showing the machine learning pipeline from data loading through Random Forest classification and biomarker export</desc>
  <defs>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="context-stroke" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  <mask id="imagine-text-gaps-01cmb0" maskUnits="userSpaceOnUse"><rect x="0" y="0" width="680" height="780" fill="white"/><rect x="101.625" y="31.25" width="166.75" height="21.5" fill="black" rx="2"/><rect x="110.6953125" y="50.5" width="148.609375" height="19" fill="black" rx="2"/><rect x="433.13671875" y="31.25" width="123.7265625" height="21.5" fill="black" rx="2"/><rect x="414.140625" y="50.5" width="161.71875" height="19" fill="black" rx="2"/><rect x="250.984375" y="151.25" width="178.03125" height="21.5" fill="black" rx="2"/><rect x="218.8359375" y="170.5" width="242.328125" height="19" fill="black" rx="2"/><rect x="269.16015625" y="235.25" width="141.6796875" height="21.5" fill="black" rx="2"/><rect x="227.19921875" y="254.5" width="225.6015625" height="19" fill="black" rx="2"/><rect x="254.53125" y="319.25" width="170.9375" height="21.5" fill="black" rx="2"/><rect x="252.11328125" y="340.5" width="175.7734375" height="19" fill="black" rx="2"/><rect x="231.87890625" y="358.5" width="216.2421875" height="19" fill="black" rx="2"/><rect x="268.0546875" y="415.25" width="143.890625" height="21.5" fill="black" rx="2"/><rect x="226.796875" y="436.5" width="226.40625" height="19" fill="black" rx="2"/><rect x="244.87890625" y="454.5" width="190.2421875" height="19" fill="black" rx="2"/><rect x="130.4375" y="511.25" width="109.125" height="21.5" fill="black" rx="2"/><rect x="101.7109375" y="532.5" width="166.578125" height="19" fill="black" rx="2"/><rect x="143.68359375" y="550.5" width="82.6328125" height="19" fill="black" rx="2"/><rect x="425.17578125" y="511.25" width="139.6484375" height="21.5" fill="black" rx="2"/><rect x="423.49609375" y="532.5" width="143.0078125" height="19" fill="black" rx="2"/><rect x="426.50390625" y="550.5" width="136.9921875" height="19" fill="black" rx="2"/><rect x="288.34375" y="643.25" width="103.3125" height="21.5" fill="black" rx="2"/><rect x="237.30078125" y="664.5" width="205.3984375" height="19" fill="black" rx="2"/><rect x="235.2265625" y="682.5" width="209.546875" height="19" fill="black" rx="2"/><rect x="26.1484375" y="39" width="19.703125" height="19" fill="black" rx="2"/><rect x="24.890625" y="159" width="22.21875" height="19" fill="black" rx="2"/><rect x="25.03515625" y="243" width="21.9296875" height="19" fill="black" rx="2"/><rect x="24.703125" y="331" width="22.59375" height="19" fill="black" rx="2"/><rect x="25" y="427" width="22" height="19" fill="black" rx="2"/><rect x="24.96484375" y="523" width="22.0703125" height="19" fill="black" rx="2"/><rect x="25.1796875" y="655" width="21.640625" height="19" fill="black" rx="2"/></mask></defs>

  <!-- Stage 1: Input data (two boxes) -->
  <g onclick="sendPrompt('What is in the normalized gene expression matrix?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="60" y="20" width="250" height="56" rx="10" stroke-width="0.5" style="fill:rgb(225, 245, 238);stroke:rgb(15, 110, 86);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="185" y="42" text-anchor="middle" dominant-baseline="central" style="fill:rgb(8, 80, 65);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Gene expression matrix</text>
    <text x="185" y="60" text-anchor="middle" dominant-baseline="central" style="fill:rgb(15, 110, 86);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">luad_normalizedData.csv</text>
  </g>

  <g onclick="sendPrompt('What does the TCGA LUAD metadata contain?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="370" y="20" width="250" height="56" rx="10" stroke-width="0.5" style="fill:rgb(225, 245, 238);stroke:rgb(15, 110, 86);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="495" y="42" text-anchor="middle" dominant-baseline="central" style="fill:rgb(8, 80, 65);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Sample metadata</text>
    <text x="495" y="60" text-anchor="middle" dominant-baseline="central" style="fill:rgb(15, 110, 86);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">TCGA_LUAD_metadata.csv</text>
  </g>

  <!-- Both converge to alignment step -->
  <path d="M185 76 L185 116 L340 116 L340 140" fill="none" stroke="#1D9E75" stroke-width="1.5" marker-end="url(#arrow)" style="fill:none;stroke:rgb(29, 158, 117);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
  <path d="M495 76 L495 116 L340 116" fill="none" stroke="#1D9E75" stroke-width="1.5" style="fill:none;stroke:rgb(29, 158, 117);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 2: Align -->
  <g onclick="sendPrompt('Why is barcode alignment important?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="160" y="140" width="360" height="56" rx="10" stroke-width="0.5" style="fill:rgb(230, 241, 251);stroke:rgb(24, 95, 165);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="340" y="162" text-anchor="middle" dominant-baseline="central" style="fill:rgb(12, 68, 124);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Align samples by barcode</text>
    <text x="340" y="180" text-anchor="middle" dominant-baseline="central" style="fill:rgb(24, 95, 165);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Intersect barcodes · subset both datasets</text>
  </g>

  <line x1="340" y1="196" x2="340" y2="224" stroke="#378ADD" marker-end="url(#arrow)" style="fill:none;stroke:rgb(115, 114, 108);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 3: Define target + split -->
  <g onclick="sendPrompt('What is the classification target in this model?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="160" y="224" width="360" height="56" rx="10" stroke-width="0.5" style="fill:rgb(230, 241, 251);stroke:rgb(24, 95, 165);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="340" y="246" text-anchor="middle" dominant-baseline="central" style="fill:rgb(12, 68, 124);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Define target &amp; split</text>
    <text x="340" y="264" text-anchor="middle" dominant-baseline="central" style="fill:rgb(24, 95, 165);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Tumor vs Normal · 80% train / 20% test</text>
  </g>

  <line x1="340" y1="280" x2="340" y2="308" stroke="#378ADD" marker-end="url(#arrow)" style="fill:none;stroke:rgb(115, 114, 108);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 4: RFE -->
  <g onclick="sendPrompt('How does Recursive Feature Elimination select genes?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="140" y="308" width="400" height="68" rx="10" stroke-width="0.5" style="fill:rgb(238, 237, 254);stroke:rgb(83, 74, 183);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="340" y="330" text-anchor="middle" dominant-baseline="central" style="fill:rgb(60, 52, 137);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Feature selection — RFE</text>
    <text x="340" y="350" text-anchor="middle" dominant-baseline="central" style="fill:rgb(83, 74, 183);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Recursive Feature Elimination</text>
    <text x="340" y="368" text-anchor="middle" dominant-baseline="central" style="fill:rgb(83, 74, 183);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Top 100 genes selected from ~20,000</text>
  </g>

  <line x1="340" y1="376" x2="340" y2="404" stroke="#7F77DD" marker-end="url(#arrow)" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(115, 114, 108);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 5: Train -->
  <g onclick="sendPrompt('What are the Random Forest hyperparameters used?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="140" y="404" width="400" height="68" rx="10" stroke-width="0.5" style="fill:rgb(238, 237, 254);stroke:rgb(83, 74, 183);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="340" y="426" text-anchor="middle" dominant-baseline="central" style="fill:rgb(60, 52, 137);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Train Random Forest</text>
    <text x="340" y="446" text-anchor="middle" dominant-baseline="central" style="fill:rgb(83, 74, 183);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">n_estimators = 100 · random_state = 42</text>
    <text x="340" y="464" text-anchor="middle" dominant-baseline="central" style="fill:rgb(83, 74, 183);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Fit on 100 RFE-selected features</text>
  </g>

  <line x1="340" y1="472" x2="340" y2="500" stroke="#7F77DD" marker-end="url(#arrow)" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(115, 114, 108);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 6: Evaluate (two outputs side by side) -->
  <g onclick="sendPrompt('What was the model accuracy on the test set?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="60" y="500" width="250" height="68" rx="10" stroke-width="0.5" style="fill:rgb(250, 236, 231);stroke:rgb(153, 60, 29);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="185" y="522" text-anchor="middle" dominant-baseline="central" style="fill:rgb(113, 43, 19);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Evaluate model</text>
    <text x="185" y="542" text-anchor="middle" dominant-baseline="central" style="fill:rgb(153, 60, 29);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Accuracy · Confusion matrix</text>
    <text x="185" y="560" text-anchor="middle" dominant-baseline="central" style="fill:rgb(153, 60, 29);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">F1-score: 1.00</text>
  </g>

  <g onclick="sendPrompt('What do the top feature importances represent?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="370" y="500" width="250" height="68" rx="10" stroke-width="0.5" style="fill:rgb(250, 236, 231);stroke:rgb(153, 60, 29);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="495" y="522" text-anchor="middle" dominant-baseline="central" style="fill:rgb(113, 43, 19);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Feature importance</text>
    <text x="495" y="542" text-anchor="middle" dominant-baseline="central" style="fill:rgb(153, 60, 29);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Gini importance ranking</text>
    <text x="495" y="560" text-anchor="middle" dominant-baseline="central" style="fill:rgb(153, 60, 29);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">Top 10 genes visualized</text>
  </g>

  <path d="M340 472 L340 486 L185 486 L185 500" fill="none" stroke="#7F77DD" stroke-width="1.5" marker-end="url(#arrow)" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(127, 119, 221);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
  <path d="M340 472 L340 486 L495 486 L495 500" fill="none" stroke="#7F77DD" stroke-width="1.5" marker-end="url(#arrow)" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(127, 119, 221);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Converge to export -->
  <path d="M185 568 L185 614 L340 614 L340 632" fill="none" stroke="#D85A30" stroke-width="1.5" marker-end="url(#arrow)" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(216, 90, 48);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
  <path d="M495 568 L495 614 L340 614" fill="none" stroke="#D85A30" stroke-width="1.5" mask="url(#imagine-text-gaps-01cmb0)" style="fill:none;stroke:rgb(216, 90, 48);color:rgb(0, 0, 0);stroke-width:1.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>

  <!-- Stage 7: Export -->
  <g onclick="sendPrompt('What files are exported from the ML pipeline?')" style="fill:rgb(0, 0, 0);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto">
    <rect x="140" y="632" width="400" height="68" rx="10" stroke-width="0.5" style="fill:rgb(250, 238, 218);stroke:rgb(133, 79, 11);color:rgb(0, 0, 0);stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:16px;font-weight:400;text-anchor:start;dominant-baseline:auto"/>
    <text x="340" y="654" text-anchor="middle" dominant-baseline="central" style="fill:rgb(99, 56, 6);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:14px;font-weight:500;text-anchor:middle;dominant-baseline:central">Export results</text>
    <text x="340" y="674" text-anchor="middle" dominant-baseline="central" style="fill:rgb(133, 79, 11);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">selected_genes.csv · ml-results.csv</text>
    <text x="340" y="692" text-anchor="middle" dominant-baseline="central" style="fill:rgb(133, 79, 11);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:central">100 biomarker candidates identified</text>
  </g>

  <!-- Step labels on left -->
  <text x="36" y="53" text-anchor="middle" fill="#1D9E75" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">01</text>
  <text x="36" y="173" text-anchor="middle" fill="#378ADD" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">02</text>
  <text x="36" y="257" text-anchor="middle" fill="#378ADD" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">03</text>
  <text x="36" y="345" text-anchor="middle" fill="#7F77DD" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">04</text>
  <text x="36" y="441" text-anchor="middle" fill="#7F77DD" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">05</text>
  <text x="36" y="537" text-anchor="middle" fill="#D85A30" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">06</text>
  <text x="36" y="669" text-anchor="middle" fill="#BA7517" style="fill:rgb(61, 61, 58);stroke:none;color:rgb(0, 0, 0);stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;opacity:1;font-family:&quot;Anthropic Sans&quot;, -apple-system, &quot;system-ui&quot;, &quot;Segoe UI&quot;, sans-serif;font-size:12px;font-weight:400;text-anchor:middle;dominant-baseline:auto">07</text>
</svg>


Raw Gene Expression Data
        │
        ▼
Metadata Alignment (Barcode Indexing)
        │
        ▼
Common Sample Subsetting
        │
        ▼
Target Label Definition
        │
        ▼
Feature Selection (Recursive Feature Elimination - RFE)
        │
        ▼
Train/Test Split
        │
        ▼
Random Forest Classifier Training
        │
        ▼
Model Evaluation (Accuracy, Confusion Matrix, F1-Score)
        │
        ▼
Biomarker Gene Identification
```

---

## ⚙️ Methods

### 1. Data Preparation

- **Barcode indexing:** Metadata barcodes are set as the DataFrame index to ensure alignment with the gene expression matrix.
- **Sample subsetting:** Only samples present in both the metadata and gene expression dataset are retained.
- **Target definition:** The classification target is binary — Primary Tumor (1) or Solid Tissue Normal (0).

### 2. Feature Selection — Recursive Feature Elimination (RFE)

RFE iteratively removes the least important features by training a model and pruning weak predictors. This step:
- Reduces the high-dimensional gene expression space to a compact set of informative genes
- Removes noise and irrelevant features that could cause overfitting
- Produces a ranked list of genes most predictive of tumor vs. normal status

### 3. Model Training — Random Forest Classifier

A Random Forest ensemble of decision trees was trained on the RFE-selected features:

| Parameter | Value |
|-----------|-------|
| Algorithm | Random Forest |
| Task | Binary Classification |
| Input | RFE-selected gene expression features |
| Output | Primary Tumor / Solid Tissue Normal |

Random Forest is well-suited for this task because it:
- Handles high-dimensional biological data effectively
- Provides feature importance scores for biomarker ranking
- Is robust to overfitting compared to single decision trees

### 4. Model Evaluation

The classifier was evaluated using:
- **Accuracy Score**
- **Confusion Matrix**
- **Precision, Recall, and F1-Score** (per class)

---

## 📊 Results

| Metric | Primary Tumor | Solid Tissue Normal |
|--------|--------------|---------------------|
| Precision | 1.00 | 1.00 |
| Recall | 1.00 | 1.00 |
| F1-Score | 1.00 | 1.00 |
| **Overall Accuracy** | **1.00** | |

The Random Forest classifier achieved **perfect classification accuracy (100%)** on the LUAD dataset with zero misclassifications.

### Confusion Matrix

```
                     Predicted
                   Tumor | Normal
Actual  Tumor  [   TP   |   0   ]
        Normal [    0   |   TN  ]
```

> ⚠️ **Important Note on Overfitting:** The perfect accuracy should be interpreted cautiously given the small dataset size (n=40). Cross-validation on a larger cohort is recommended to validate generalizability.

---

## 🧬 Biological Significance

Genes selected by RFE as top features represent candidate biomarkers that:
- Distinguish tumor from normal lung tissue at the transcriptomic level
- May be involved in oncogenic pathways or tumor suppression
- Are consistent with findings from the parallel DEG analysis (upregulated/downregulated gene sets)

These features, when cross-referenced with differentially expressed genes, can help narrow the list of clinically relevant biomarker candidates.

---

## 📁 Output Files

| File | Description |
|------|-------------|
| `ml-results.csv` | Classification results per sample |
| `ML_Lung_Cancer.ipynb` | Full analysis notebook |

---

## 🛠️ Requirements

```bash
pip install pandas numpy scikit-learn matplotlib seaborn jupyter
```

| Package | Purpose |
|---------|---------|
| `pandas` | Data manipulation and indexing |
| `numpy` | Numerical operations |
| `scikit-learn` | RFE, Random Forest, evaluation metrics |
| `matplotlib` / `seaborn` | Visualization |

---

## ▶️ How to Run

```bash
# Clone the repository
git clone https://github.com/sanzidaanee/Hackbio-cancer-internship.git

# Navigate to the Stage 3 code directory
cd "Hackbio-cancer-internship/Stage 3/Code"

# Launch Jupyter Notebook
jupyter notebook ML_Lung_Cancer.ipynb
```

---

## 📚 References

- Abdelwahab O, Awad N, Elserafy M, Badr E (2022). A feature selection-based framework to identify biomarkers for cancer diagnosis: A focus on lung adenocarcinoma. *PLOS ONE* 17(9): e0269126.
- The Cancer Genome Atlas (2015). TCGA Data Portal. National Cancer Institute.

---

## 👥 Contributors

Sanzida Akhter Anee, Sk Arif, Nada Ghozlan, Mennatallah Mohamed Ebrahim Mahmoud, Stéphie Raveloson, Chidimma Nwaku, Zaka Ullah, Idahosa Clinton

**HackBio Cancer Bioinformatics Internship — Stage 3**

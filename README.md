# Gittributor / Attributor

A mass spectrometry expert tool developed at **IPAG (Institut de Planétologie et d'Astrophysique de Grenoble, France)**.

## 📥 Quick Download

**Direct download links:**
- [ATTRIBUTOR 2022.pxp](ATTRIBUTOR%202022.pxp) (~55 MB) - Main project file
- [MainProc.ipf](MainProc.ipf) (~813 KB) - Source code
- [Igor Pro Official Website](https://www.wavemetrics.com/) - WaveMetrics

## 📚 Documentation

A comprehensive user guide is available in French:
- **[Attributor - Complete User Guide (French)](attributor-guide-full.html)**

This guide covers:
- Installation and launch procedures
- Detailed user interface description
- Step-by-step workflow
- Advanced features and algorithms
- Mathematical foundations
- Practical tutorials
- Troubleshooting
- License information

## 🚀 Quick Start

### Prerequisites
- **Igor Pro 64-bit** (WaveMetrics) must be installed
- Compatible with **Windows 7/10/11**

### Installation
1. Ensure both files are in the same directory:
   - `ATTRIBUTOR 2022.pxp` (main project file)
   - `MainProc.ipf` (source code)

2. Launch **Igor64.exe**

3. Open the project: **File → Open → Open Experiment...** → Select `ATTRIBUTOR 2022.pxp`

### Basic Usage
- Use the **periodic table panel** to click on elements and build molecular formulas
- The isotopic distribution simulation appears automatically in the **elaborateur** window
- Compare with experimental data
- Use **Recal** button to calibrate experimental data
- Add simulations to **agregador** for multi-comparison

## 📁 Files

| File | Description |
|------|-------------|
| `ATTRIBUTOR 2022.pxp` | Main Igor Pro project file |
| `MainProc.ipf` | Source code with all functions and procedures |
| `attributor-guide-full.html` | Complete user guide (French, all-in-one HTML) |
| `index.html` | Redirect page |

## 🔬 Features

- **Interactive molecular formula building** via periodic table
- **Isotopic distribution simulation** with high precision
- **Comparison** of simulated spectra with experimental data
- **Calibration** of experimental data
- **Multi-molecule comparison** via agregador
- **Mass defect analysis** via DMVM window
- **MS/MS mode** support

## 📊 Main Windows

- **panel**: Periodic table interface for building molecules
- **elaborateur**: Main graph window for visualization
- **agregador**: Comparison window for multiple simulations
- **dmvm**: Mass defect analysis window
- **AdvancedManager**: Data management and loading

## 🎓 Tutorial Example

To identify an unknown molecule with a peak at m/z ≈ 214.08485:
1. Click 10× C, 10× H, 4× O in the periodic table
2. Formula C10H10O4 appears
3. Isotopic simulation is displayed
4. Compare with experimental data
5. Adjust if needed and validate

## 🐛 Troubleshooting

- **"Function not found"**: Ensure MainProc.ipf is in the same directory
- **Windows not opening**: Check Igor Pro History window for errors
- **Slow simulation**: Increase isotopic combo probability threshold
- **Mass mismatch**: Use Recal button, check ion charge

## 📜 License

Non-exclusive, non-transferable license for research purposes only.

© 2014-2026 François-Régis ORTHOUS-DAUNAY, IPAG (UMR5374), CNRS, UJF Grenoble, France

## 🏢 About IPAG

**Institut de Planétologie et d'Astrophysique de Grenoble**
- UMR 5274 (CNRS / Université Grenoble Alpes)
- Website: [https://ipag.osug.fr](https://ipag.osug.fr)
- Research domains: Astrophysics, Planetology, Earth and Universe Sciences, Instrumentation

## 📞 Contact

Developer: François-Régis ORTHOUS-DAUNAY

---

*Last updated: September 2026*

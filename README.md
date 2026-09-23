# Attributor / Gittributor

An expert mass spectrometry analysis tool developed at **IPAG (Institut de Planetologie et d'Astrophysique de Grenoble, France)**.

## Quick Download

**Direct download links:**
- [ATTRIBUTOR 2022.pxp](ATTRIBUTOR%202022.pxp) (~55 MB) - Main project file
- [MainProc.ipf](MainProc.ipf) (~813 KB) - Source code
- [Igor Pro Official Website](https://www.wavemetrics.com/) - WaveMetrics

## Documentation

A comprehensive expert user guide is available:
- **[Attributor - Complete Expert User Guide](attributor-guide-full.html)**

This guide covers:
- Scientific foundations and theory
- System requirements and prerequisites
- Installation and launch procedures
- Detailed user interface description
- Standard and advanced workflows
- Mathematical algorithms and calculations
- Step-by-step tutorials
- Technical reference and Igor Pro commands
- Troubleshooting guide
- License information

## Quick Start

### Prerequisites
- **Igor Pro 64-bit** (WaveMetrics) version 8.x or higher - <strong>32-bit NOT supported</strong>
- Compatible with **Windows 7/10/11** (64-bit)
- Minimum **4 GB RAM** recommended, **8-16 GB** for large molecules

### Installation
1. Ensure both files are in the same directory:
   - `ATTRIBUTOR 2022.pxp` (main project file)
   - `MainProc.ipf` (source code with all procedures)

2. Launch **Igor64.exe** (64-bit version only)

3. Open the project: **File -> Open -> Open Experiment...** -> Select `ATTRIBUTOR 2022.pxp`

   **Alternative:** Drag and drop the .pxp file onto Igor64.exe

### Basic Usage
- Use the **periodic table panel** to click on elements and build molecular formulas
- The isotopic distribution simulation appears automatically in the **elaborateur** window
- Red lines = Theoretical simulation, Gray lines = Experimental data
- Compare with experimental data using **Recal** button for automatic calibration
- Add simulations to **agregador** for multi-molecule comparison

## Files

| File | Description | Size | Required |
|------|-------------|------|----------|
| `ATTRIBUTOR 2022.pxp` | Main Igor Pro experiment file | ~55 MB | Yes |
| `MainProc.ipf` | Source code with all functions and procedures | ~813 KB | Yes |
| `attributor-guide-full.html` | Complete expert user guide (all-in-one HTML) | ~150 KB | No |
| `index.html` | Redirect page | ~1 KB | No |

## Features

### Core Capabilities
- **Interactive molecular formula building** via periodic table interface
- **High-precision isotopic distribution simulation** using polynomial generator method
- **Real-time experimental comparison** with logarithmic intensity scaling
- **Mass calibration** (automatic and manual modes)
- **Multi-molecule comparison** via aggregator window
- **Mass defect analysis** (DMVM plots, Kendrick analysis)
- **MS/MS mode** support for tandem mass spectrometry

### Advanced Features
- Fast Fourier Transform (FFT) optimized calculations
- Support for all stable isotopes with natural abundances
- Dynamic range: 0.001% to 100% of base peak
- Multiple peak profile shapes (Gaussian, Lorentzian, hybrid)
- Region Of Interest (ROI) selection
- Batch processing capabilities

### Window System
- **panel**: Periodic table control panel
- **elaborateur**: Main graph window for visualization
- **agregador**: Multi-simulation comparison window
- **dmvm**: Mass defect vs m/z analysis window
- **AdvancedManager**: Data management and spectrum loading

## Tutorial Example

To identify an unknown molecule with a peak at m/z = 214.08485:
1. Use AdvancedManager to load experimental spectrum
2. In dmvm window, note mass defect: +0.08485 Da
3. Build formula C10H10O4 in panel (10x C, 10x H, 4x O)
4. Monoisotopic mass: 214.05790880 Da
5. Compare isotopic profile in elaborateur
6. Use Recal if peaks do not align perfectly
7. Add to agregador to compare with alternative hypotheses

## Troubleshooting

| Problem | Cause | Solution |
|--------|-------|----------|
| Igor Pro fails to launch | Corrupted installation | Reinstall Igor Pro 64-bit |
| Project fails to load | MainProc.ipf not in same directory | Place both files together |
| "Function not found" error | Procedure not compiled | Verify MainProc.ipf location |
| Simulation is slow | Large molecule (>100 atoms) | Increase threshold or reduce size |
| Peaks do not align | Mass calibration needed | Use Recal button |
| Memory errors | Insufficient RAM | Close other applications or increase RAM |

## Technical Details

### System Requirements
- **CPU:** Dual-core 2 GHz (minimum), Quad-core 3 GHz (recommended), Multi-core >3.5 GHz (optimal)
- **RAM:** 4 GB (minimum), 8 GB (recommended), 16+ GB (optimal for large molecules)
- **Storage:** 100 MB free space
- **Graphics:** OpenGL-compatible GPU recommended
- **Display:** 1280x720 (minimum), 1920x1080 (recommended)

### Igor Pro Compatibility
| Version | Status | Notes |
|--------|--------|-------|
| Igor Pro 9.x | Fully Supported | Optimized for new features |
| Igor Pro 8.x | Fully Supported | All features available |
| Igor Pro 7.x (64-bit) | Limited Support | Some UI issues possible |
| Igor Pro 6.x or lower | Not Supported | Incompatible API |
| Igor Pro 32-bit | Not Supported | Memory limitations |

## License

Non-exclusive, non-transferable license for **research purposes only**.

### Permitted Uses
- Academic research in mass spectrometry and related fields
- Non-commercial scientific analysis
- Educational purposes in accredited institutions
- Collaborative research with IPAG or affiliated groups

### Restricted Uses
- Commercial applications without prior written consent
- Redistribution of software or source code
- Modification of source code without permission
- Use in for-profit organizations without appropriate licensing

### Copyright

&copy; 2014-2026 **Francois-Regis ORTHOUS-DAUNAY**

**Institution:** IPAG (Institut de Planetologie et d'Astrophysique de Grenoble)

**Affiliations:** UMR 5274, CNRS, Universite Grenoble Alpes, France

## About IPAG

**Institut de Planetologie et d'Astrophysique de Grenoble**

- **UMR 5274** (CNRS / Universite Grenoble Alpes)
- **Website:** [https://ipag.osug.fr](https://ipag.osug.fr)
- **Research domains:** Astrophysics, Planetology, Earth and Universe Sciences, Instrumentation

## Contact

Developer: **Francois-Regis ORTHOUS-DAUNAY**

---

*Last updated: September 2026*

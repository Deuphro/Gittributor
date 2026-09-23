# ATTRIBUTOR - Mass Spectrometry Formula Attribution Tool

[![WaveMetrics Igor Pro](https://www.wavemetrics.com/images/wavemetrics_logo.png)](https://www.wavemetrics.com/)

ATTRIBUTOR is an advanced software tool for **molecular formula attribution** from high-resolution mass spectrometry data. Developed at **IPAG (Institut de Planétologie et d'Astrophysique de Grenoble)** by François-Régis ORTHOUS-DAUNAY, this Igor Pro-based application enables researchers to analyze isotopic patterns and mass-to-charge (m/z) ratios to identify and validate molecular compositions.

## Features

### Molecular Formula Construction
- **Interactive element selection**: Click buttons for each element to build molecular formulas
- **Comprehensive element support**: All periodic table elements available
- **Charge state management**: Handle positive and negative ion states
- **Visual formula representation**: Real-time display of molecular composition

### Isotopic Pattern Analysis
- **Multi-isotope combination calculations**: Computes all possible isotopic combinations
- **Probability-weighted distributions**: Theoretical isotopic abundance calculations
- **Mass defect analysis**: Precise mass-to-charge ratio calculations
- **Isotope filtering**: Adjustable thresholds for probability criteria

### Spectral Analysis Tools
- **Data visualization**: Multi-panel graph displays with customizable parameters
- **Peak detection**: Automatic identification of spectral peaks
- **Noise filtering**: Advanced algorithms for signal-to-noise enhancement
- **Pattern matching**: Compare theoretical and experimental spectra

### Advanced Capabilities
- **Graph customization**: Logarithmic and linear scaling options
- **Profile fitting**: Gaussian, Lorentzian, and asymmetric peak profiles
- **Data export**: Save and export analysis results
- **Multi-spectrum overlay**: Compare multiple datasets simultaneously

## Installation

### Requirements
- **Igor Pro 8.00 or later** (required for long wave name support)
- Windows operating system
- Minimum 4 GB RAM (8 GB recommended for large datasets)

### Quick Start

1. **Download the files**:
   - [ATTRIBUTOR 2022.pxp](ATTRIBUTOR%202022.pxp) - Main experiment file
   - [MainProc.ipf](MainProc.ipf) - Core procedure file

2. **Place files in the same directory**

3. **Launch ATTRIBUTOR**:
   ```bash
   Igor64.exe ".\ATTRIBUTOR 2022.pxp"
   ```

## Usage

### Basic Workflow

1. **Load your mass spectrometry data** into Igor Pro
2. **Open ATTRIBUTOR** using the command above
3. **Build your molecular formula** by clicking element buttons
4. **Adjust isotopic criteria** using the threshold controls
5. **Compare theoretical patterns** with your experimental data
6. **Validate and refine** your molecular attribution

### Element Selection

The ATTRIBUTOR interface provides buttons for all elements organized by periodic table groups:
- **Light elements**: H, He, Li, Be, B, C, N, O, F, Ne
- **Alkali/Alkaline Earth**: Na, Mg, Al, Si, P, S, Cl, Ar
- **Transition metals**: Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, etc.
- **Heavy elements**: All remaining elements up to U

### Isotopic Simulation

ATTRIBUTOR calculates isotopic distributions using:
- Element-specific isotopic masses and abundances
- Binomial probability distributions
- Configurable probability thresholds (0-100%)
- Charge state corrections

### Graph Windows

ATTRIBUTOR uses multiple graph windows:
- **elaborateur**: Main analysis window with data overlay
- **agregateur**: Aggregated spectrum display
- **dmvm**: Mass defect visualization
- **molmanager**: Molecular formula management

## File Descriptions

### ATTRIBUTOR 2022.pxp
- **Type**: Igor Pro Experiment file (packed)
- **Size**: ~57.3 MB
- **Content**: All experiment data, windows, variables, and saved analysis
- **Note**: Binary format, requires Igor Pro to open

### MainProc.ipf
- **Type**: Igor Pro Procedure file
- **Size**: ~831 KB
- **Content**: Core functions and macros including:
  - `genestringmol()` - Generate molecular formula strings
  - `choisylecrible()` - Select isotope filtering method
  - `generesimu()` - Generate isotopic simulations
  - `crible1-10()` - Isotope combination calculations
  - Element button handlers (elemH, elemC, elemN, etc.)
  - Graph window management functions

## Technical Details

### Supported Elements
ATTRIBUTOR supports all naturally occurring elements with their isotopic compositions:
- Up to 10 isotopes per element
- Accurate isotopic masses and natural abundances
- Configurable for custom isotopic data

### Mass Range
- Theoretical mass range: 0-10,000 m/z (configurable)
- Practical range depends on instrument capabilities
- High-resolution mass defect calculations

### Algorithm
- **Combinatorial approach**: Exhaustive calculation of all possible isotopic combinations
- **Probability weighting**: Based on natural isotopic abundances
- **Threshold filtering**: Configurable minimum probability for inclusion
- **Charge correction**: Electron mass subtraction/addition

## License

ATTRIBUTOR is provided under a **non-exclusive, non-transferable license** for research activities only:

- **Permitted**: Use for personal research activities, archiving, and backup
- **Restricted**: 
  - Distribution to third parties
  - Sub-licensing or copying
  - Reverse engineering or modification
  - Commercial use or service provision
  - Third-party access

**All rights reserved**. The software, documentation, and specifications remain the intellectual property of the owner (François-Régis ORTHOUS-DAUNAY).

## Author

**François-Régis ORTHOUS-DAUNAY**

- **Affiliation**: IPAG (UMR5374), CNRS, Université Joseph Fourier
- **Location**: Grenoble, France
- **Date**: September 1st, 2014 (original version)

## References

- **IPAG**: [Institut de Planétologie et d'Astrophysique de Grenoble](https://www.ipag.fr/)
- **CNRS**: [Centre National de la Recherche Scientifique](https://www.cnrs.fr/)
- **Igor Pro**: [WaveMetrics Igor Pro](https://www.wavemetrics.com/)

## Version History

- **2022**: Current version with updated features
- **2014**: Original release (September 1st)

## Support

For questions or issues regarding ATTRIBUTOR:

1. Ensure you have **Igor Pro 8.00 or later** installed
2. Verify both files (.pxp and .ipf) are in the same directory
3. Check that you have sufficient memory for your dataset
4. Consult the WaveMetrics Igor Pro documentation

## Web Interface

A companion web page is available with:
- Complete feature overview
- Direct download links
- Usage instructions
- Official links to WaveMetrics

Open [index.html](index.html) in your browser for the full documentation.

---

*ATTRIBUTOR is provided "as is" without warranty of any kind. Use at your own risk.*

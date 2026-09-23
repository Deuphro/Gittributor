# ATTRIBUTOR

**ATTRIBUTOR** is an Igor Pro research application for exploring molecular-formula hypotheses from high-resolution mass-spectrometry data. It combines bounded formula search, charge-aware isotope simulation, formula-list management, graph-based propagation, mass-defect visualization, and Hough-style pattern discovery.

This repository contains a compact English companion website and the source files included with the application.

## Files

| File | Size | Role |
|---|---:|---|
| [`ATTRIBUTOR 2022.pxp`](ATTRIBUTOR%202022.pxp) | 57,343,327 bytes (57.3 MB) | Packed Igor Pro experiment: the primary file to open. |
| [`MainProc.ipf`](MainProc.ipf) | 831,617 bytes (831 KB) | UTF-8 Igor procedure source accompanying the experiment. |
| [`index.html`](index.html) | — | Self-contained English feature and usage page. |

## Requirements

- Windows
- A licensed [Igor Pro](https://www.wavemetrics.com/products/igorpro/) installation capable of opening the supplied experiment
- The Igor Pro XOPs referenced by the source (`SaveGraph`, `Resize Controls`, and `Multi-peak fitting 2.0`) when compiling the `.ipf` directly

The repository does not document a minimum Igor Pro version. Use a current Igor Pro release and consult WaveMetrics' platform-support documentation for current operating-system requirements.

## Launch

From PowerShell, in the directory containing the downloaded files:

```powershell
Igor64.exe '.\ATTRIBUTOR 2022.pxp'
```

If `Igor64.exe` is not available on `PATH`, replace it with the full path to the executable in your Igor Pro installation.

Keep the `.pxp` and `.ipf` together if you also plan to inspect or compile the procedure source.

## Typical workflow

1. **Load a spectrum** from an X–Y mass/intensity text file.
2. **Activate it** in `AdvancedManager`, or load it into the current ROI.
3. **Build a formula** with the `mendeleiev` element panel, formula fields, and ion-charge control.
4. **Inspect its isotope envelope** in `elaborateur`, where sticks, an asymmetric profile, mass labels, and ppm labels are overlaid on the spectrum.
5. **Add useful patterns to a factor list** and set minimum/maximum stoichiometric bounds in `attributeur`.
6. **Calculate or batch-attribute** target peaks, then review the proposals through their complete isotope simulations.
7. **Compare and save formula lists** in `molmanager`; overlay them in `agregateur` and inspect them in `dmvm`.
8. **Explore recurring families** with Graph Attribution and the Hough transform when a compositional pattern is not known in advance.

A formula produced by ATTRIBUTOR is a candidate hypothesis, not an identification by itself.

## Main windows

- **`mendeleiev`** — builds stoichiometries, sets the ion charge, controls the relative isotope-combination threshold, and sends formulas to factors or the aggregator.
- **`elaborateur`** — compares measured and simulated isotope envelopes; exposes profile width, shape, asymmetry, FWHM, and Gaussian/Lorentz shortcuts.
- **`agregateur`** — overlays experimental spectra and multiple formula simulations; supports log scaling, sticks/profile mode, labels, focus, rescaling, and formula-related context actions.
- **`attributeur`** — sets a target m/z, selects stored factors and count constraints, and launches single, batch, or graph-based attribution.
- **`AdvancedManager`** — catalogs spectra and chromatograms; controls activation, copies, removal, and spectrum/mass-defect overlays.
- **`molmanager`** — stores reusable formula bags and exposes mass, ppm bias, focus, rescaling, ROI conversion, and list operations.
- **`dmvm` / `datastat`** — plots mass defect against m/z and uses cumulative intensity/point distributions to build ROIs or launch Hough analysis.
- **`filtres`** — provides Union-Find peak detection, radio-noise rejection, anti-ringing, noise cuts, normalization, convolution, calibration, and ROI utilities.
- **`PeakQualCont`** — manages peak profiles, summed fits, residuals, baselines, and reusable peak sets.
- **`AdvMolManPanel` / `AdvMolManGraph` / `AdvMolManOper`** — synchronize formula-list selection, configurable projections, clustering, ratio plots, and K-naries views.
- **`etudeTransformationHough`** — displays the slope/intercept Hough accumulator, marginal histograms, reference-composition lines, zoom tools, and region-to-ROI/attribution actions.

The operation window includes a **Partitioning** tab, but the supplied build does not implement controls for it.

## Computational methods

### Bounded formula search — `trouvemass`

The target m/z is converted to a neutral-mass target using the selected charge. Stored factors are represented by their most abundant-isotope masses, then sorted. A cumulative-product mixed-radix walk enumerates bounded stoichiometric vectors. All but the final, heaviest factor are enumerated explicitly; that final factor is solved from the remaining mass error. Branches known to overshoot are skipped.

The output is ranked by absolute signed ppm error. Constraint bounds, charge, and the selected factor collection define the search space; the algorithm is not an unrestricted elemental search.

### Isotope simulation — `genestringmol`

For each element, ATTRIBUTOR enumerates isotope occupancies across up to ten modeled isotope states. Their probability follows the multinomial form:

```text
P(k₀ … kᵣ) = n! / ∏ kᵢ! · ∏ pᵢᵏⁱ
Σ kᵢ = n
```

Local mass/probability pairs are combined across elements. A relative probability threshold prunes weak local combinations. The complete envelope is charge-corrected, scaled to the active spectral context, and used to generate stick and profile views plus mass-defect traces.

### Graph Attribution — `kruskal4mass`

Graph Attribution compares measured peak pairs with mass gaps taken from selected reference compositions. Valid pairs become weighted edges. `kruskal4mass` processes edges in increasing distance order and uses union–find to build connected mass families while enforcing a maximum node degree. A formula can then be seeded at the current target or a family root and propagated depth-first by adding or removing each edge's stoichiometry and charge.

### DMVM and the Hough transform

The mass-defect coordinate is:

```text
d(m) = m − round(m)
```

For two DMVM points, ATTRIBUTOR derives the line parameters:

```text
slope     = (dⱼ − dᵢ) / (mⱼ − mᵢ)
intercept = dᵢ − slope · mᵢ
```

Every pair votes in a slope/intercept accumulator. Bright cells and marginal histograms reveal linear structures supported by many pairs. Reference formula lines can be overlaid when known, but the first discovery is data-driven: it can expose a recurring family before its chemistry is named. A selected region can then seed an ROI or attribution workflow.

## License

The header of `MainProc.ipf` grants a non-exclusive, non-transferable license for the recipient's research use of the software and documentation, and for an archival/backup copy of ATTRIBUTOR, provided that titles, trademarks, and other personal or reserved rights are retained on copies and remain subject to the license.

The recipient must not make the software usable by or distribute it to third parties, sublicense, copy, reverse engineer, dissociate/decompose, or modify it, or use it for service exchange, subcontracting, a service environment, or third-party access. Copyright, trade secrets, and other intellectual-property rights remain with the owner.

Read the complete header in [`MainProc.ipf`](MainProc.ipf) before reuse or redistribution.

## Attribution

Developed by **François-Régis ORTHOUS-DAUNAY** at **IPAG (UMR5374), CNRS, UJF Grenoble, France**. The source header is dated September 2014; the supplied packed experiment is named `ATTRIBUTOR 2022.pxp`.

## Web page

Open [`index.html`](index.html) directly in a browser, or serve the repository with any static web server. The page has no external JavaScript, font, image, or stylesheet dependencies.

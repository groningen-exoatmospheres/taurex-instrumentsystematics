# ☁️ TauREx-InstrumentSystematics: A TauREx plugin to handle instrumental systematics.

A custom [TauREx](https://github.com/ucl-exoplanets/TauREx3_public) plugin to correct for the main instrumental systematics: includes offsets, slopes, and instrument response functions. 
---
## 💾 Installation 

You can install TauREx-InstrumentSystematics using python in the plugin's folder:

```python
python3 setup.py install
```
---

## 📌 Purpose

This module can be used to remove instrument systematics, particularly when using multiple instruments. It can also be used to load instrument files containing the instrument's reponse function to be convolved with the theoretical TauREx model before fit.

---

## 🔧 Model Parameters

The plugin is used by specifying "spectra_instr" as keyword in the parfile.

| Name | Description |
|------|-------------|
| `path_spectra` | Path to spectra. These should be .txt with format: wave (um), transit depth, error, wave width (um)  |
| `offsets` | Array of same size as path_spectra with offsets |
| `slopes` | Array of same size as path_spectra with slopes |
| `broadening_profile` | Path to instrument files. These should be .fits using STScI formal or .txt using formal: wavelength (um), resolution |
| `wlshift` | To introduce a linear wavelength shift. THIS THIS NOT YET IMPLEMENTED AND SO DOES NOTHING. FOR WAVELENGTH SHIFTS USE TAUREX-MULTIMODEL FOR NOW. |  

---





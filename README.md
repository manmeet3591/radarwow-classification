# radarwow-classification

**radar-wow stage 5: Hydrometeor classification and melting layer.** Melting-layer detection from QVPs, fuzzy-logic / clustering classes.

Part of [radar-wow](https://github.com/radar-wow), a set of repos that together build a radar QPE system that works for S, C and X band radars.
Every algorithm is one class implementing one paper, with the same interface: it takes an xradar sweep or volume with canonical moment names, declares the bands and moments it needs, and records provenance.

## Sub-components

| sub-component | algorithm | source | DOI | status |
|---|---|---|---|---|
| quasi-vertical profiles | `quasi_vertical_profile` | Ryzhkov et al. 2016, JTECH | [doi:10.1175/JTECH-D-15-0020.1](https://doi.org/10.1175/JTECH-D-15-0020.1) | implemented + notebook |
| melting layer from QVPs | `detect_melting_layer` | Sanchez-Rivas & Rico-Ramirez 2021, AMT | [doi:10.5194/amt-14-2873-2021](https://doi.org/10.5194/amt-14-2873-2021) | implemented + notebook |
| temporal consistency | `MeltingLayerTracker` | radar-wow: rejects jumps > 0.5 km, optional model prior |  | implemented + notebook |
| NEXRAD MLDA | `Giangrande2008ML` | Giangrande, Krause & Ryzhkov 2008, JAMC | [doi:10.1175/2007JAMC1634.1](https://doi.org/10.1175/2007JAMC1634.1) | planned |
| fuzzy-logic HCA | `Park2009HCA` | Park et al. 2009, Wea. Forecasting | [doi:10.1175/2008WAF2222205.1](https://doi.org/10.1175/2008WAF2222205.1) | planned |
| semi-supervised clustering | `Besic2016` | Besic et al. 2016, AMT | [doi:10.5194/amt-9-4425-2016](https://doi.org/10.5194/amt-9-4425-2016) | planned |
| non-meteorological echo (dust, bio) | `DualPolNonMetClassifier` | Kilambi et al. 2018; Gourley et al. 2007 (see survey/01) |  | planned |

## Install

```bash
pip install git+https://github.com/radar-wow/radarwow-core git+https://github.com/radar-wow/radarwow-classification
```

## Notebooks (free data, S / C / X band)

- The melting layer is used and checked in the [radarwow-verify KTLX notebook](https://github.com/radar-wow/radarwow-verify/blob/main/notebooks/01_ktlx_20220505_vs_mrms.ipynb)

## Usage

```python
from radarwow.core import load_sample
from radarwow.classification import *  # see table above

vol = load_sample('arm_csapr2_c')        # or open_volume('my_file')
# out = SomeAlgorithm()(vol)            # band read from the file, or band='S'|'C'|'X'
```

## Status legend

*implemented + notebook*: code, tests, executed notebook. *planned*: interface fixed, code to write. *stub*: class exists, raises `NotImplementedError`. *to review*: paper to read before deciding.

## License

MIT. Data used in the notebooks comes from public archives (NEXRAD on AWS, ARM, IMD via open-radar-data) under their own terms.

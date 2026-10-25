# radarwow-classification

**radar-wow stage 5: Hydrometeor classification and melting layer.** Fuzzy-logic / clustering classes, melting-layer detection.

Part of [radar-wow](https://github.com/radar-wow), a set of repos that together build a radar QPE system that works for S, C and X band radars.
Every algorithm is one class implementing one paper, with the same interface: it takes an xradar sweep or volume with canonical moment names, declares the bands and moments it needs, and records provenance.

## Sub-components

| sub-component | algorithm | source | DOI | status |
|---|---|---|---|---|
| fuzzy-logic HCA | `Park2009HCA` | Park et al. 2009, Wea. Forecasting | [doi:10.1175/2008WAF2222205.1](https://doi.org/10.1175/2008WAF2222205.1) | planned |
| melting layer | `Giangrande2008ML` | Giangrande, Krause & Ryzhkov 2008, JAMC | [doi:10.1175/2007JAMC1634.1](https://doi.org/10.1175/2007JAMC1634.1) | planned |
| semi-supervised clustering | `Besic2016` | Besic et al. 2016, AMT | [doi:10.5194/amt-9-4425-2016](https://doi.org/10.5194/amt-9-4425-2016) | planned |
| biological echo class | see `radarwow-qc` `MistNet` | Lin et al. 2019 | [doi:10.1111/2041-210X.13280](https://doi.org/10.1111/2041-210X.13280) | stub |

## Install

```bash
pip install git+https://github.com/radar-wow/radarwow-core git+https://github.com/radar-wow/radarwow-classification
```

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

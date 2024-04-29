# SKUF
Compressor for image. Now support BMP image format and numpy SVD ("npSVD"),
power method SVD ("pwmSVD"), block power method SVD ("bpmSVD") compress algorithms.
Compress image to [SKUF](SKUF-description.md).

---
Эксперимент [тутъ](experiment.md)

---

## Usage
### CLI
Implemented CLI for SKUF lib

To see the usage, run:
```Bash
python src/compressor.py --help
```

### SKUFlib

Main class for SKUFlib is SKUF. Init it with image or skuf file.
If you init with image, it automatically compress it by given algo.

Example:
```python
from src.SKUFlib.SKUF import *
from src.SKUFlib.SVD import *

skuf = SKUF("img_src/10-3-1.bmp", BlockPowerMethodSVD(compression_degree=3, tolerance=10e-2))
skuf.save("img_src/10-3-1.skuf")
skuf_load = SKUF("img_src/10-3-1.skuf")
img = skuf_load.decode()
img.save("img_src/10-3-1_decompress.skuf")
```


# Small Knowledge Unified Format!

A SKUF file is intended for storing compressed RGB 24 bit images. At the moment only SVD algorithms are implemented,
but the format supports other algorithms correspond to the Transformer interface.

File description:

| Name                      | Format | Size in bytes   | Value                                                   |
|---------------------------|--------|-----------------|---------------------------------------------------------|
| magic                     | <4s    | 4               | b'skuf'                                                 |
| version                   | <b     | 1               | 1                                                       |
| size x                    | <I     | 4               | width of original image in pixels                       |
| size y                    | <I     | 4               | height of original image in pixels                      |
| format                    | <8s    | 8               | original image format (only BMP supported now)          |
| algo                      | <8s    | 8               | algo used for compression ("npSVD", "pwmSVD", "bpmSVD") |
| color data size           | <I     | 4               | size of algo depending data per color                   |
| algo depending data red   | -      | color data size | description for SVD below                               |
| algo depending data green | -      | color data size | description for SVD below                               |
| algo depending data blue  | -      | color data size | description for SVD below                               |

SVD data:

| Name        | Format | Size in bytes      | Value                     |
|-------------|--------|--------------------|---------------------------|
| height      | <I     | 4                  | height of original image  |
| saving part | <I     | 4                  | saving part of SVD matrix |
| width       | <I     | 4                  | width of original image   |
| u           | -      | height*saving_part | u matrix                  |
| s           | -      | saving_part        | s vector                  |
| v           | -      | width*saving_part  | v matrix                  |


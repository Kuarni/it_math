from copy import deepcopy
from struct import pack, unpack, calcsize
from typing import NamedTuple, BinaryIO

import numpy as np
from PIL import Image

from .Transformer import Transformer
from .SVD import NpSVD, PowerMethodSVD, BlockPowerMethodSVD
from .Packing import Packing


class SKUF:
    """Small Knowledge Unified Format"""

    class __ImgMetadata(NamedTuple):
        size_x: int
        size_y: int
        format: bytes
        algo: bytes
        size_of_one_color_data: int

    __str_byte_sizes = 8
    __magic: Packing = Packing(b'skuf', "<4s")
    __version: Packing = Packing(1, "<b")
    __img_metadata: Packing = Packing(None, f"<II{__str_byte_sizes}s{__str_byte_sizes}sI")
    supported_formats = ["BMP"]
    supported_algos = {
        NpSVD.name: NpSVD,
        PowerMethodSVD.name: PowerMethodSVD,
        BlockPowerMethodSVD.name: BlockPowerMethodSVD
    }

    def __init__(self, file_path: str, encoder: Transformer = None):
        """Load scuf or compress given image to scuf via encoder"""
        self.__supported_formats_b = map(self.__str_to_bytes, self.supported_formats)
        self.__supported_algos_b = {self.__str_to_bytes(k): v for k, v in self.supported_algos.items()}
        self.__skuf_r = None
        self.__skuf_g = None
        self.__skuf_b = None
        with open(file_path, "rb") as f:
            if unpack(self.__magic.packing, f.read(calcsize(self.__magic.packing)))[0] == self.__magic.value:
                self.load(f)
            else:
                self.encode(Image.open(file_path), encoder)

    def __str_to_bytes(self, s: str) -> bytes:
        b = bytes(s, encoding='ascii')
        if len(b) > self.__str_byte_sizes:
            b = b[:self.__str_byte_sizes]
        else:
            b = b + b'\0' * (self.__str_byte_sizes - len(b))
        return b

    def encode(self, img, encoder: Transformer):
        if encoder is None:
            encoder = self.supported_algos["npSVD"]()
        if encoder.name not in self.supported_algos:
            raise ValueError(
                f"Unsupported encoder type {encoder.name}\nList of supported encoders: {self.supported_algos.keys()}")
        if img.format not in self.supported_formats:
            raise ValueError(
                f"Unsupported image format \"{img.format}\"\nList of supported formats: {self.supported_formats}")

        r, g, b = img.split()
        self.__skuf_r = deepcopy(encoder).encode(r)
        self.__skuf_g = deepcopy(encoder).encode(g)
        self.__skuf_b = deepcopy(encoder).encode(b)
        self.__img_metadata.value = self.__ImgMetadata(img.size[0], img.size[1], self.__str_to_bytes(img.format),
                                                       self.__str_to_bytes(encoder.name),
                                                       self.__skuf_r.size())

    def decode(self):
        if self.__skuf_r is None or self.__skuf_g is None or self.__skuf_b is None or self.__img_metadata.value is None:
            raise ValueError("Image is empty, nothing to decode")
        img = np.zeros((self.__img_metadata.value.size_y, self.__img_metadata.value.size_x, 3), dtype=np.uint8)
        img[:, :, 0] = self.__skuf_r.decode()
        img[:, :, 1] = self.__skuf_g.decode()
        img[:, :, 2] = self.__skuf_b.decode()
        return Image.fromarray(img)

    def load(self, file: BinaryIO, offset=0):
        file.seek(offset)
        name = unpack(self.__magic.packing, file.read(calcsize(self.__magic.packing)))[0]
        if name != self.__magic.value:
            raise ValueError("Not a skuf file was given")
        version = unpack(self.__version.packing, file.read(calcsize(self.__version.packing)))[0]
        if version != self.__version.value:
            raise ValueError("Unsupported version")
        self.__img_metadata.value = self.__ImgMetadata(
            *unpack(self.__img_metadata.packing, file.read(calcsize(self.__img_metadata.packing))))
        if self.__img_metadata.value.format not in self.__supported_formats_b:
            raise ValueError(f"Unsupported image format \"{self.__img_metadata.value.format}\"")
        if self.__img_metadata.value.algo not in self.__supported_algos_b:
            raise ValueError(
                f"Unsupported image compression algorithm \"{self.__img_metadata.value.algo}\" in skuf file")
        algo = self.__supported_algos_b[self.__img_metadata.value.algo]
        self.__skuf_r = algo().read(file)
        self.__skuf_g = algo().read(file)
        self.__skuf_b = algo().read(file)

    def save(self, file_path: str):
        with open(file_path, 'wb') as file:
            file.write(pack(self.__magic.packing, self.__magic.value))
            file.write(pack(self.__version.packing, self.__version.value))
            file.write(pack(self.__img_metadata.packing, *self.__img_metadata.value))
            file.write(self.__skuf_r.to_bytes())
            file.write(self.__skuf_g.to_bytes())
            file.write(self.__skuf_b.to_bytes())

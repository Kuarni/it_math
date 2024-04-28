from abc import abstractmethod
from typing import BinaryIO

import PIL.Image


class Transformer:
    name: str

    @abstractmethod
    def encode(self, image: PIL.Image):
        pass

    @abstractmethod
    def size(self):
        pass

    @abstractmethod
    def decode(self):
        pass

    @abstractmethod
    def to_bytes(self):
        pass

    @abstractmethod
    def read(self, f: BinaryIO):
        pass

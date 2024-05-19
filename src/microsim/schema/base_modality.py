from typing import Literal
from typing import Protocol

import numpy as np

from microsim.schema._base_model import SimBaseModel

MODALITIES = Literal["widefield", "confocal", "two-photon", "light-sheet"]


class ModalityProtocol(Protocol):
    type: MODALITIES
    def generate_psf(self, *args, **kwargs) -> np.typing.ArrayLike: ...
    def render(self, truth, channel, objective_lens, settings, xp) -> None: ...

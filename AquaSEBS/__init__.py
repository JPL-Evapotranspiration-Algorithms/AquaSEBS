from .AquaSEBS import AquaSEBS
from .water_heat_flux import water_heat_flux
from .exceptions import MissingOfflineParameter
from .constants import RESAMPLING_METHOD, MASK_NON_WATER_PIXELS, UPSCALE_TO_DAYLIGHT
from .version import __version__

__author__ = "Gregory H. Halverson"

__all__ = [
    "AquaSEBS",
    "water_heat_flux",
    "MissingOfflineParameter",
    "RESAMPLING_METHOD",
    "MASK_NON_WATER_PIXELS",
    "UPSCALE_TO_DAYLIGHT",
    "__version__",
    "__author__",
    "__citation__",
    "__bibtex__",
]

__citation__ = (
    f"Halverson, G. H. (2026). AquaSEBS: Aquatic Surface Energy Balance System "
    f"Water Surface Evaporation Model (v{__version__}). "
    f"https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS"
)

__bibtex__ = f"""@software{{halverson2026aquasebs,
  title={{AquaSEBS: Aquatic Surface Energy Balance System Water Surface Evaporation Model}},
  author={{Halverson, Gregory H.}},
  version={{{__version__}}},
  year={{2026}},
  publisher={{GitHub}},
  url={{https://github.com/JPL-Evapotranspiration-Algorithms/AquaSEBS}}
}}

@article{{abdelrady2016aquasebs,
  title={{Surface Energy Balance of Fresh and Saline Waters: AquaSEBS}},
  author={{Abdelrady, Ahmed and Timmermans, Joris and Vekerdy, Zolt{{\'a}}n and Salama, Mhd. Suhyb}},
  journal={{Remote Sensing}},
  volume={{8}},
  number={{7}},
  pages={{583}},
  year={{2016}},
  publisher={{MDPI}},
  doi={{10.3390/rs8070583}}
}}

@article{{fisher2023remotely,
  title={{Remotely sensed terrestrial open water evaporation}},
  author={{Fisher, Joshua B and Dohlen, Matthew B and Halverson, Gregory H and Collison, Jacob W and Pearson, Christopher and Huntington, Justin L}},
  journal={{Scientific Reports}},
  volume={{13}},
  number={{1}},
  pages={{8174}},
  year={{2023}},
  publisher={{Nature Publishing Group}},
  doi={{10.1038/s41598-023-34921-2}}
}}"""

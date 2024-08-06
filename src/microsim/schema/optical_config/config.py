import inspect
from collections.abc import Sequence
from typing import Any, cast

import xarray as xr
from pydantic import Field, model_validator
from scipy.constants import c, h

from microsim.fpbase import SpectrumOwner
from microsim.schema._base_model import SimBaseModel
from microsim.schema.sample.fluorophore import Fluorophore
from microsim.schema.spectrum import Spectrum

from .filter import Filter, Placement, SpectrumFilter


class LightSource(SimBaseModel):
    name: str = ""
    spectrum: Spectrum
    power: float | None = None  # W/cm^2

    @classmethod
    def from_fpbase(cls, light: SpectrumOwner) -> "LightSource":
        return cls(name=light.name, spectrum=Spectrum.from_fpbase(light.spectrum))

    def plot(self, show: bool = True) -> None:
        self.spectrum.plot(show=show)

    @classmethod
    def laser(cls, wavelength: float, power: float | None = None) -> "LightSource":
        return cls(
            name=f"{wavelength}nm Laser",
            spectrum=Spectrum(wavelength=[wavelength], intensity=[1]),
            power=power,
        )


class OpticalConfig(SimBaseModel):
    name: str = ""
    filters: list[Filter] = Field(default_factory=list)
    lights: list[LightSource] = Field(default_factory=list)
    # seemingly duplicate of power in LightSource
    # but it depends on where the power is being measured
    # TODO: it's tough deciding where power should go...
    # it could also go on Simulation itself as a function of space.
    power: float | None = None  # total power of all lights after filters

    def absorption_rate(self, fluorophore: Fluorophore) -> "xr.DataArray":
        """Return the absorption rate of a fluorophore in this configuration.

        The absorption rate is the number of photons absorbed per second per
        fluorophore.
        """
        # get irradiance scaled to power
        irrad = self.irradiance  # W/cm^2
        # absorption cross section in cm^2
        cross_section = fluorophore.absorption_cross_section
        # calculate excitation rate (this takes care of finding overlapping wavelengths)
        watts_absorbed = irrad * cross_section
        wavelength_meters = cast("xr.DataArray", watts_absorbed.coords["w"] * 1e-9)
        joules_per_photon = h * c / wavelength_meters
        abs_rate = watts_absorbed / joules_per_photon  # 1/s

        # add metadata
        abs_rate.name = "absorption_rate"
        abs_rate.attrs["long_name"] = "Absorption rate"
        abs_rate.attrs["units"] = "photons/s/fluorophore"
        return abs_rate  # type: ignore [no-any-return]

    @property
    def excitation(self) -> Filter | None:
        """Combine all excitation filters into a single spectrum."""
        filters = []
        for f in self.filters:
            if f.placement in {Placement.EX_PATH, Placement.BS_INV, Placement.ALL}:
                filters.append(f)
            if f.placement == Placement.BS:
                filters.append(f.inverted())
        return self._merge(filters, spectrum="excitation")

    @property
    def illumination(self) -> Spectrum | None:
        exc = self.excitation
        if self.lights:
            l0, *rest = self.lights
            illum_spect = l0.spectrum
            if rest:
                for light in rest:
                    illum_spect = illum_spect * light.spectrum
            if exc:
                return illum_spect * exc.spectrum
            return illum_spect
        return exc.spectrum if exc else None

    @property
    def irradiance(self) -> "xr.DataArray":
        if (illum := self.illumination) is None:
            raise ValueError("This Optical Config has no illumination spectrum.")
        # get irradiance scaled to power
        irrad = illum.as_xarray()  # W/cm^2
        # normalize area under curve to 1
        irrad = irrad / irrad.sum()
        # scale to power
        # if self.power is not None:
        # irrad = irrad * self.power
        irrad.name = "irradiance"
        irrad.attrs["long_name"] = "Irradiance"
        irrad.attrs["units"] = "W/cm^2"
        return irrad

    @property
    def emission(self) -> Filter | None:
        """Combine all emission filters into a single spectrum."""
        filters = []
        for f in self.filters:
            if f.placement in {Placement.EM_PATH, Placement.BS, Placement.ALL}:
                filters.append(f)
            if f.placement == Placement.BS_INV:
                filters.append(f.inverted())
        return self._merge(filters, spectrum="emission")

    def _merge(
        self, filters: Sequence[Filter], spectrum: str = "spectrum"
    ) -> Filter | None:
        if not filters:
            return None
        if len(filters) == 1:
            return filters[0]
        effective_spectrum = filters[0].spectrum
        for filt in filters[1:]:
            effective_spectrum = effective_spectrum * filt.spectrum
        return SpectrumFilter(
            name=f"Effective {spectrum} for {self.name}",
            transmission=effective_spectrum,
        )

    @classmethod
    def from_fpbase(
        cls, microscope_id: str, config_name: str | None = None
    ) -> "OpticalConfig":
        from microsim.fpbase import get_microscope

        if config_name is None:
            if "::" not in microscope_id:  # pragma: no cover
                raise ValueError(
                    "If config_name is not provided, microscope_id must be "
                    "in the form 'scope::config'"
                )
            microscope_id, config_name = microscope_id.split("::")

        fpbase_scope = get_microscope(microscope_id)
        for cfg in fpbase_scope.opticalConfigs:
            if cfg.name.lower() == config_name.lower():
                if cfg.light:
                    lights = [LightSource.from_fpbase(cfg.light)]
                else:
                    lights = []
                return cls(
                    name=cfg.name,
                    filters=[SpectrumFilter.from_fpbase(f) for f in cfg.filters],
                    lights=lights,
                )

        raise ValueError(
            f"Could not find config named {config_name!r} in FPbase microscope "
            f"{microscope_id!r}. Available names: "
            f"{', '.join(repr(c.name) for c in fpbase_scope.opticalConfigs)}"
        )

    @model_validator(mode="before")
    def _vmodel(cls, value: Any) -> Any:
        if isinstance(value, str):
            if "::" not in value:  # pragma: no cover
                raise ValueError(
                    "If OpticalConfig is provided as a string, it must be "
                    "in the form 'fpbase_scope_id::config_name'"
                )
            # TODO: seems weird to have to cast back to dict...
            # but otherwise doesn't work with 'before' validator.  look into it.
            return cls.from_fpbase(value).model_dump()
        return value

    def plot(self, show: bool = True) -> None:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(12, 3))
        ax = fig.add_subplot(111)

        legend = []
        for filt in self.filters:
            ax.plot(filt.spectrum.wavelength, filt.spectrum.intensity)
            legend.append(filt.name)
        if any(legend):
            ax.legend(legend)
        if show:
            plt.show()

    def all_spectra(self) -> "xr.DataArray":
        data, coords = [], []
        for filt in self.filters:
            data.append(filt.spectrum.as_xarray())
            coords.append(f"{filt.name} ({filt.placement.name})")
        for light in self.lights:
            data.append(light.spectrum.as_xarray())
            coords.append(light.name)
        da: xr.DataArray = xr.concat(data, dim="spectra")
        da.coords.update({"spectra": coords})
        return da

    # WARNING: dark magic ahead
    # This is a hack to make OpticalConfig hashable and comparable, but only
    # when used in the context of a pandas DataFrame or xarray DataArray coordinate.
    # this allows syntax like `data_array.sel(c='FITC')` to work as expected.
    def __hash__(self) -> int:
        frame = inspect.stack()[1]
        if "pandas" in frame.filename and frame.function == "get_loc":
            return hash(self.name)
        return id(self)

    def __eq__(self, value: object) -> bool:
        frame = inspect.stack()[1]
        if "pandas" in frame.filename and frame.function == "get_loc":
            return hash(self.name) == hash(value)
        return super().__eq__(value)

    def __str__(self):
        return self.name

"""Read and write geodesy_modeling Insar2dObject instances as HDF5.

HDF5 format v1 (one object per file):

    /lon, /lat                        1D coordinate datasets
    /LOS, /LOS_unc                    2D observation/uncertainty datasets
    /lkv_E, /lkv_N, /lkv_U            2D ground-to-satellite look vectors
    /coherence                       2D coherence dataset

All 2D arrays have shape (len(lat), len(lon)), with southern rows first.
The root attributes contain format_name, format_version, look_direction,
and optional starttime/endtime (ISO 8601 strings), plus optional los_units.

Requires numpy and h5py. Place this file alongside class_model.py in
geodesy_modeling/datatypes/InSAR_2D_Object, or import it as a separate module.
"""

from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

from geodesy_modeling.datatypes.InSAR_2D_Object.class_model import Insar2dObject


FORMAT_NAME = "geodesy_modeling.Insar2dObject"
FORMAT_VERSION = 1
ARRAY_FIELDS = ("lon", "lat", "LOS", "LOS_unc", "lkv_E", "lkv_N", "lkv_U", "coherence")
GRID_FIELDS = ("LOS", "LOS_unc", "lkv_E", "lkv_N", "lkv_U", "coherence")


def _check_arrays(arrays):
    """Check the schema, including dimensions the class constructor does not check."""
    for name in ARRAY_FIELDS:
        arr = arrays[name]
        if arr.dtype.kind not in "biufc":
            raise TypeError(f"{name} must contain numeric data, not {arr.dtype}")
        if arr.ndim != (1 if name in ("lon", "lat") else 2):
            raise ValueError(f"{name} has the wrong number of dimensions: {arr.shape}")

    if len(arrays["lat"]) < 2 or len(arrays["lon"]) < 1:
        raise ValueError("The grid needs at least two latitudes and one longitude")

    expected = (len(arrays["lat"]), len(arrays["lon"]))
    for name in GRID_FIELDS:
        if arrays[name].shape != expected:
            raise ValueError(f"{name} has shape {arrays[name].shape}, expected {expected}")


def _write_optional_datetime(h5file, name, value):
    if value is None:
        return  # Missing attribute means None.
    if not isinstance(value, datetime):
        raise TypeError(f"{name} must be a datetime.datetime or None")
    h5file.attrs[name] = value.isoformat()


def _read_optional_datetime(h5file, name):
    value = h5file.attrs.get(name)
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    try:
        return datetime.fromisoformat(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {name} timestamp in HDF5 file: {value!r}") from exc


def write_insar2d_hdf5(insar, filename, *, overwrite=False, compression=4, los_units=None,
                       verbose=False):
    """Save an Insar2dObject to a single, lossless, self-describing HDF5 file.

    Parameters
    ----------
    insar : Insar2dObject
        Object to write.
    filename : str or pathlib.Path
        Destination .h5 file.
    overwrite : bool, default False
        Refuse to replace an existing file unless explicitly enabled.
    compression : int or None, default 4
        Gzip level 1--9 for 2D arrays; 0/None disables compression.
        Compression is lossless; NaNs and array dtypes are preserved.
    los_units : str or None, default None
        Optional metadata (e.g., "mm" or "radians"). The class does not
        itself store units, so no unit is assumed here.
    verbose : True or False, default False
        Print the file path or not.

    Returns
    -------
    pathlib.Path
        Path of the saved file.
    """
    if verbose:
        print("Writing file %s" % filename)
    if not isinstance(insar, Insar2dObject):
        raise TypeError("insar must be an Insar2dObject")
    if compression is not None and (not isinstance(compression, int) or not 0 <= compression <= 9):
        raise ValueError("compression must be an integer from 0 to 9 or None")
    if los_units is not None and not isinstance(los_units, str):
        raise TypeError("los_units must be a string or None")
    if insar.look_direction not in ("right", "left"):
        raise ValueError("look_direction must be 'right' or 'left'")
    for name in ("starttime", "endtime"):
        value = getattr(insar, name)
        if value is not None and not isinstance(value, datetime):
            raise TypeError(f"{name} must be a datetime.datetime or None")

    arrays = {}
    for name in ARRAY_FIELDS:
        field = getattr(insar, name)
        if np.ma.isMaskedArray(field) and np.ma.getmaskarray(field).any():
            raise ValueError(f"{name} is masked; fill masked cells (e.g., with NaN) before writing")
        arrays[name] = np.asarray(field)
    _check_arrays(arrays)

    path = Path(filename)
    mode = "w" if overwrite else "x"
    with h5py.File(path, mode) as h5file:
        h5file.attrs["format_name"] = FORMAT_NAME
        h5file.attrs["format_version"] = FORMAT_VERSION
        h5file.attrs["look_direction"] = insar.look_direction
        h5file.attrs["grid_order"] = "lat_lon"  # rows south-to-north, columns west-to-east
        h5file.attrs["look_vector_direction"] = "ground_to_satellite"
        _write_optional_datetime(h5file, "starttime", insar.starttime)
        _write_optional_datetime(h5file, "endtime", insar.endtime)
        if los_units is not None:
            h5file.attrs["los_units"] = los_units

        for name in ARRAY_FIELDS:
            kwargs = {}
            if name in GRID_FIELDS and compression:
                kwargs = {"compression": "gzip", "compression_opts": compression,
                          "shuffle": True, "fletcher32": True}
            dataset = h5file.create_dataset(name, data=arrays[name], **kwargs)
            if name == "lon":
                dataset.attrs["units"] = "degrees_east"
            elif name == "lat":
                dataset.attrs["units"] = "degrees_north"
            elif name in ("lkv_E", "lkv_N", "lkv_U", "coherence"):
                dataset.attrs["units"] = "1"
            elif los_units is not None:
                dataset.attrs["units"] = los_units
    return path


def read_insar2d_hdf5(filename, verbose=False):
    """Read an HDF5 format-v1 file and reconstruct an Insar2dObject.

    Parameters
    ----------
    filename : string
        Path to the hdf5 filename
    verbose : True or False, default False
        Print the file path or not.

    Returns
    -------
    Insar2dObject
    """
    if verbose:
        print("Reading file %s" % filename)
    with h5py.File(filename, "r") as h5file:
        format_name = h5file.attrs.get("format_name")
        if format_name != FORMAT_NAME:
            raise ValueError(f"Not an {FORMAT_NAME} HDF5 file: {format_name!r}")
        version = h5file.attrs.get("format_version")
        if version != FORMAT_VERSION:
            raise ValueError(f"Unsupported Insar2dObject HDF5 format version: {version!r}")

        arrays = {}
        for name in ARRAY_FIELDS:
            if name not in h5file or not isinstance(h5file[name], h5py.Dataset):
                raise ValueError(f"Missing HDF5 dataset: /{name}")
            arrays[name] = h5file[name][...]
        _check_arrays(arrays)

        look_direction = h5file.attrs.get("look_direction")
        if isinstance(look_direction, bytes):
            look_direction = look_direction.decode("utf-8")
        starttime = _read_optional_datetime(h5file, "starttime")
        endtime = _read_optional_datetime(h5file, "endtime")

    # Construct after closing the file: all arrays are ordinary NumPy arrays.
    return Insar2dObject(**arrays, starttime=starttime, endtime=endtime,
                         look_direction=look_direction)

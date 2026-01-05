#!/usr/bin/env python3
"""
brukerread.py
==============

Simple Python replacement for the MATLAB function `brukerread.m`
that reads Bruker NMR data files (.dsc/.dta and .par/.spc).

Only the standard library and numpy are required.
"""

import os
import sys
import numpy as np

# --------------------------------------------------------------------------- #
# Utility helpers – mimic a few MATLAB conveniences
# --------------------------------------------------------------------------- #

def safeget(dct, key, default):
    """Return dictionary[key] if present, otherwise default."""
    return dct.get(key, default)

def isfield(dct, key):
    return key in dct

# --------------------------------------------------------------------------- #
# Parsing the .dsc / .par text files
# --------------------------------------------------------------------------- #

def parse_dsc(filename):
    """
    Read a Bruker .dsc or .par file.
    Each line has the form  KEY : VALUE.
    Comments (lines starting with # or ! ) and empty lines are ignored.
    Returns a dictionary of key/value pairs (both strings).
    """
    if not os.path.isfile(filename):
        raise IOError(f"File not found: {file_path}")

    # 2ï¸â£  Load the entire file in one go
    with open(filename, "r", encoding="utf-8", errors="ignore") as fh:
        raw_data = fh.read()
    
    # 3ï¸â£  Split into lines (Pythonâs builtâin str.splitlines)
    lines = raw_data.splitlines()
    
    # 4ï¸â£  Parse
    params = {}
    previous_key = None
    
    forbidden = '~!@#$%^&*().//\\'
    
    for raw_line in lines:
        line = raw_line.strip()
        if not line:                 # skip empty lines
            continue
    
        # Skip lines that start with a forbidden character
        if line[0] in forbidden:
            continue
            
        # Split at the first whitespace
        parts = line.split(None, 1)
        if not parts:                # malformed line â ignore
            continue
    
        key, rest = parts[0], parts[1] if len(parts) > 1 else ""
    
        # New field â key does NOT start with a digit
        if not key[0].isdigit():
            previous_key = key
            params[key] = rest.strip().replace("'", '') 
        # Continuation line â key starts with a digit
        else:
            if previous_key is None:  # no preceding key â ignore
                continue
            existing = params.get(previous_key, "")
            params[previous_key] = (existing + rest).strip()
    
    return params

    # try:
    #     with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
    #         for line in f:
    #             line = line.strip()
    #             if not line or line.startswith('#') or line.startswith('!'):
    #                 continue
    #             if ':' in line:
    #                 key, val = line.split(':', 1)
    #                 dsc[key.strip()] = val.strip()
    # except FileNotFoundError:
    #     # In the MATLAB code the function will simply return an empty dsc
    #     pass
    # return dsc


def parse_par(par):
    """
    Convert a MATLAB âstructâ (dictionary) of parameters into the
    dictionary expected by the rest of the application.

    Parameters
    ----------
    par : dict
        Input parameters (keys are strings).

    Returns
    -------
    res : dict
        Result dictionary with axes, titles, etc.
    """
    if not isinstance(par, dict):
        raise ValueError("Input must be a dictionary of string parameters")

    res = {}

    # ----- topâlevel parameters ---------------------------------------
    if 'MWFQ' in par:
        res['freq1'] = _to_number(par['MWFQ'])
    if 'MF' in par:
        res['freq1'] = _to_number(par['MF'])
    if 'CenterField' in par:
        res['cf'] = _to_number(strip_unit(par['CenterField'])) * 1e-4

    if 'IKKF' in par:
        res['complex'] = _is_equal(par['IKKF'], 'CPLX')
    elif 'XQAD' in par:
        res['complex'] = _is_equal(par['XQAD'], 'ON')
    else:
        res['complex'] = 0

    if 'TITL' in par:
        res['title'] = par['TITL']
    elif 'JCO' in par:
        res['title'] = par['JCO']
    else:
        res['title'] = '?'

    # ----- axes -------------------------------------------------------
    res = getrespar321(par, res, 'x')
    res = getrespar321(par, res, 'y')

    return res

def _to_number(value):
    """Return a float from a string or number; raise ValueError if impossible."""
    try:
        return float(value)
    except Exception:
        raise ValueError(f"Cannot convert {value!r} to a number")


def _is_equal(a, b):
    """Caseâinsensitive string comparison."""
    return str(a).strip().upper() == b.upper()
def getrespar321(par, res, axletter):
    """
    Helper that processes one axis (x or y) and updates the result
    dictionary.

    Parameters
    ----------
    par   : dict
        Input parameters.
    res   : dict
        Result dictionary to update.
    axletter : str
        Axis letter, usually 'x' or 'y'.

    Returns
    -------
    res : dict
        Updated result dictionary.
    """
    if not isinstance(par, dict):
        raise ValueError("Input must be a dictionary of string parameters")

    # initialise placeholders
    ax_lower = axletter.lower()
    res[ax_lower] = np.array([], dtype=float)
    res[f'{ax_lower}label'] = ''

    # defaults
    dim = 1
    sf = 0.0
    step = 1.0
    label = '?'
    unit = safeget(par, f'{axletter.upper()}UNI', '?')
    unit = unit.replace("'", '')

    # XSophe style handling
    if ('JSS' not in par) and ('RES' in par or 'RRES' in par) and ('HCF' in par):
        par['JSS'] = '2'
        par['RES'] = safeget(par, 'RES', safeget(par, 'RRES', '1024'))

    # Main decision tree based on JSS value
    if 'JSS' in par:
        jss_val = int(_to_number(par['JSS']))

        # --- 1D CW sweep files ----------------------------------------
        if jss_val == 2:
            if axletter.upper() != 'Y':
                dim = int(_to_number(safeget(par, 'RES', '1024')))
                if 'GST' in par and 'GSI' in par:
                    sf = _to_number(par['GST'])
                    step = _to_number(par['GSI']) / (dim - 1)
                elif 'HCF' in par and 'HSW' in par:
                    center = _to_number(par['HCF'])
                    width = _to_number(par['HSW'])
                    step = width / (dim - 1)
                    sf = center - width / 2
                elif 'HCF' in par and 'GST' in par:
                    center = _to_number(par['HCF'])
                    width = 2 * abs(center - _to_number(par['GST']))
                    step = width / (dim - 1)
                    sf = center - width / 2
                elif 'HCF' in par and 'GSI' in par:
                    par['DOS'] = '1'
                    center = _to_number(par['HCF'])
                    width = _to_number(par['GSI'])
                    step = width / (dim - 1)
                    sf = center - width / 2
                else:
                    width = abs(center - _to_number(par['GST'])) * 2
                    step = width / (dim - 1)
                    sf = center - width / 2
                label = 'Magnetic Field, G'

        # --- ESP 380 -----------------------------------------------
        elif jss_val == 32:
            if axletter.upper() != 'Y':
                if f'{axletter.upper()}QNT' in par:
                    dim = int(_to_number(par[f'{axletter.upper()}PLS']))
                    idx = par[f'{axletter.upper()}QNT'].replace(' ', '')
                    if idx == 'Time':
                        step = 0.0
                        val = [float(v) for v in par.get('Psd5', '').split()]
                        for v in val[68:75]:          # 69:75 in MATLAB
                            if v > 0:
                                res['step'] = v
                                break
                    elif idx == 'Magn.Field':
                        if 'HCF' in par:
                            cf = _to_number(par['HCF'])
                            wd = _to_number(par['HSW']) if 'HSW' in par else abs(cf - _to_number(par['GST'])) * 2
                        else:
                            cf = _to_number(par['GST'])
                            wd = _to_number(par['GSI'])
                            cf += wd / 2
                        step = wd / (dim - 1)
                        sf = cf - wd / 2
                        label = 'Magnetic Field, G'
                    elif idx in ('RF1', 'RF2'):
                        sf = float(par.get(f'{idx}StartFreq', '0'))
                        wd = float(par.get(f'{idx}SweepWidth', '0'))
                        step = wd / (dim - 1)
                        label = f'{idx}, MHz'
                    elif idx == '1.RFSource':
                        sf = _to_number(par.get('ESF', '0'))
                        wd = _to_number(par.get('ESW', '0'))
                        step = wd / (dim - 1)
                        label = 'RF, MHz'
                else:
                    dim = int(_to_number(par.get(f'{axletter.upper()}PLS', '1')))
                    if 'JUN' in par:
                        unit = par['JUN']
                    if 'GST' in par:
                        sf = _to_number(par['GST'])
                        wd = _to_number(par['GSI'])
                        step = wd / (dim - 1)
                        label = f'?,{unit}'
                    else:
                        if 'HCF' in par:
                            cf = _to_number(par['HCF']) * 1e-4
                        dim = int(_to_number(par.get(f'{axletter.upper()}PLS', '1')))
                        val = [float(v) for v in par.get('XPD9', '').split()]
                        step = val[5] * 8
                        label = 'Time, ns'

        # --- 2D files -----------------------------------------------
        elif jss_val in (4128, 4144):
            dim = int(_to_number(par.get(f'SS{axletter.upper()}', '1')))
            if axletter.upper() == 'X' and res.get('complex', 0):
                dim //= 2
            sw = _to_number(safeget(par, f'X{axletter.upper()}WI', str(dim)))
            step = sw / (dim - 1)
            unit = safeget(par, f'X{axletter.upper()}UN', '?')
            unit = unit.replace("'", '')
            label = f'Time, {unit}'

        # --- generic (1D or 2D Bruker style) ------------------------
        else:
            if f'{axletter.upper()}TYP' in par:
                if par[f'{axletter.upper()}TYP'] != 'NODATA':
                    nam = safeget(par, f'{axletter.upper()}NAM', '?')
                    dim = int(_to_number(par[f'{axletter.upper()}PTS']))
                    sf = _to_number(par[f'{axletter.upper()}MIN'])
                    wd = _to_number(par[f'{axletter.upper()}WID'])
                    step = wd / (dim - 1)
                    label = f'{nam} {unit}'
            elif f'{axletter.upper()}NAM' in par:
                dim = int(_to_number(par[f'{axletter.upper()}PTS']))
                sf = _to_number(par[f'{axletter.upper()}MIN'])
                wd = _to_number(par[f'{axletter.upper()}WID'])
                step = wd / (dim - 1)
                nam = safeget(par, f'{axletter.upper()}NAM', '?')
                label = f'{nam}, {unit}'

    # --- JEX â ESP580 -------------------------------------------------
    elif 'JEX' in par and axletter.upper() != 'Y':
        dim = int(_to_number(safeget(par, 'RES', '1024')))
        if par['JEX'] == 'ENDOR':
            sf = _to_number(par['ESF'])
            width = _to_number(par['ESW'])
            step = width / (dim - 1)
            label = 'Frequency, MHz'

    # --- Fallback when TYP not NODATA -------------------------------
    elif f'{axletter.upper()}TYP' in par and par[f'{axletter.upper()}TYP'] != 'NODATA':
        dim = int(_to_number(par[f'{axletter.upper()}PTS']))
        sf = _to_number(par[f'{axletter.upper()}MIN'])
        wd = _to_number(par[f'{axletter.upper()}WID'])
        step = wd / (dim - 1)
        nam = safeget(par, f'{axletter.upper()}NAM', '?')
        label = f'{nam}, {unit}'

        if f'{axletter.upper()}AxisQuant' in par:
            idx = par[f'{axletter.upper()}AxisQuant'].replace(' ', '')
            label = f'{idx}, {unit}'

    # ----- build the full axis array ----------------------------------
    res[ax_lower] = sf + step * np.arange(dim, dtype=float)
    res[f'{ax_lower}label'] = label

    return res

def strip_unit(s):
    """Keep only digits, a decimal point and the letter 'E' (or 'e')."""
    return ''.join(ch for ch in s if ch.isdigit() or ch in {'.', 'E', 'e'})

# --------------------------------------------------------------------------- #
# Reading the binary data matrix
# --------------------------------------------------------------------------- #

def getmatrix(filename, dims, dimorder, fmt, endian, complex_flag):
    """
    Low‑level routine that mimics MATLAB's `fread`/`reshape`/`ipermute`.

    Parameters
    ----------
    filename   : str
        Path to the binary data file.
    dims       : list of int
        Desired shape of the output array.
    dimorder   : list of int
        Order in which the dimensions are read from the file.
        Zero‑based indexing is used.
    fmt        : str
        String indicating the data type ('int32', 'float64', 'float32', 'float')
    endian     : str
        '<' for little‑endian, '>' for big‑endian.
    complex_flag : bool
        If True, the file contains interleaved real/imag parts.

    Returns
    -------
    out : numpy.ndarray
        The data array with shape `dims`.
    """
    dtype = np.dtype(endian + fmt)

    # Number of scalars to read
    N = (2 if complex_flag else 1) * np.prod(dims)

    try:
        raw = np.fromfile(filename, dtype=dtype, count=N)
    except OSError as exc:
        raise RuntimeError(f'Could not read file "{filename}": {exc}') from exc

    if raw.size < N:
        raise RuntimeError(f'Only {raw.size} elements read; {N} expected.')

    # Combine interleaved real/imag parts if needed
    if complex_flag:
        raw = raw[::2] + 1j * raw[1::2]

    # The file is stored in column‑major order with the dimensions
    # specified by `dimorder`.  We therefore first reshape using that
    # order (with Fortran ordering) and then transpose back to the
    # requested shape.
    shape_fortran = tuple(dims[i] for i in dimorder)
    out = raw.reshape(shape_fortran, order='F')
    out = out.transpose(tuple(dimorder))

    return out


# --------------------------------------------------------------------------- #
# Main routine – the Python equivalent of `brukerread`
# --------------------------------------------------------------------------- #

def brukerread(filename, return_ax=False, return_dsc=False):
    """
    Read a Bruker dataset.

    Parameters
    ----------
    filename : str
        Path to a Bruker file (extension can be .dsc/.dta or .par/.spc).
    return_ax : bool, optional
        If True, return the axis dictionary as the first output.
    return_dsc : bool, optional
        If True, return the parsed description dictionary as the last output.

    Returns
    -------
    Depending on the flags, a tuple (ax, data) or (ax, data, dsc) or just data.
    """
    # ----- Split file name and determine file pair ---------------------------------
    ppath, ext = os.path.splitext(filename)
    ext = ext.lower()

    if ext in ('.dsc', '.dta'):
        fname = ppath+'.dta'
        dscname = ppath+'.dsc'
    elif ext in ('.par', '.spc'):
        fname = ppath+'.spc'
        dscname = ppath+'.par'
    else:
        raise ValueError(f'Unsupported file extension: {ext}')

    # ----- Parse the description file ---------------------------------------------
    dsc = parse_dsc(dscname)

    # ----- Build the initial axis dictionary ---------------------------------------
    ax = parse_par(dsc)

    # ----- Default dimensions -------------------------------------------------------
    dims = [ax.get('x', np.array([], dtype=float)).size,
            ax.get('y', np.array([], dtype=float)).size, 1]

    # ----- Endianness & format ----------------------------------------------------
    endian = '<'                 # little‑endian is default
    fmt = 'f8'                   # default float64
    complex_flag = ax.get('complex', False)

    if ext in ('.dsc', '.dta'):
        if safeget(dsc, 'BSEQ', 'BIG') == 'BIG':
            endian = '>'
        irfmt = safeget(dsc, 'IRFMT', 'D')
        if irfmt == 'I':
            fmt = 'i4'
        else:
            fmt = 'f8'
        data = getmatrix(fname, dims, [0, 1, 2], fmt, endian, complex_flag)

    elif ext in ('.par', '.spc'):
        jss = int(float(safeget(dsc, 'JSS', '2')))
        if jss > 0:
            dos = int(float(safeget(dsc, 'DOS', '0')))
            if dos:
                endian = '>'
            if complex_flag:
                dims[0] *= 2
                tmp = getmatrix(fname, dims, [0, 1, 2], 'i4', endian, False)
                sz = ax.get('x', np.array([], dtype=float)).size
                data = tmp[:sz, :] + 1j * tmp[sz:, :]
            else:
                vers = int(float(safeget(dsc, 'VERS', '0')))
                fmt = 'f4' if vers != 769 else 'f8'
                data = getmatrix(fname, dims, [0, 1, 2], fmt, endian, complex_flag)
        else:
            data = np.array([], dtype=complex)
    else:
        data = np.array([], dtype=complex)

    # ----- Read optional GF axis files --------------------------------------------
    for k, typ in enumerate(('X', 'Y')):
        typ_key = f'{typ}TYP'
        if isfield(dsc, typ_key):
            axtype = dsc[typ_key]
            if axtype == 'IGD':
                gf_file = os.path.join(ppath, f'{name}.{typ}GF')
                try:
                    with open(gf_file, 'rb') as f:
                        # The GF files contain 64‑bit floats (float64)
                        ax[typ.lower()] = np.fromfile(f, dtype=np.dtype(endian + 'f8'))
                except Exception:
                    print(f'Error reading axis file: {gf_file}', file=sys.stderr)

    # ----- Fallback axes if they are still empty -------------------------------
    if ax.get('x', np.array([], dtype=float)).size != data.shape[0]:
        ax['x'] = np.arange(1, data.shape[0] + 1, dtype=float)
    if ax.get('y', np.array([], dtype=float)).size != data.shape[1]:
        ax['y'] = np.arange(1, data.shape[1] + 1, dtype=float)

    # ----- Return values ---------------------------------------------------------
    if return_ax and return_dsc:
        return ax, data, dsc
    elif return_ax:
        return ax, data
    else:
        return data


# --------------------------------------------------------------------------- #
# Example usage
# --------------------------------------------------------------------------- #

if __name__ == '__main__':
    fname = "D:\\PSU\\PSUDrive\\OneDrive - The Pennsylvania State University\\Papers\\EPR_CbHydA1\\SelectedData\\Hox\\2021_11_03\\q20211103_02.DTA" 
    # Set the flags that you want – e.g. `return_ax=True` will give you the axes
    ax, data = brukerread(fname, return_ax=True)

    print('Data shape:', data.shape)
    print('First axis (x) length:', ax['x'].size)
    print('Second axis (y) length:', ax['y'].size)

    # If you need the full description dictionary as well:
    # ax, data, dsc = brukerread(fname, return_ax=True, return_dsc=True)
    
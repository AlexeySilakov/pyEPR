#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Matlab â Python conversion of the XEPRpar routine

Only NumPy is imported; the rest uses plain Python data structures.
"""

import numpy as np


# ----------------------------------------------------------------------
# Helper functions â pure Python (no external imports)
# ----------------------------------------------------------------------
def _to_number(value):
    """Return a float from a string or number; raise ValueError if impossible."""
    try:
        return float(value)
    except Exception:
        raise ValueError(f"Cannot convert {value!r} to a number")


def _is_equal(a, b):
    """Caseâinsensitive string comparison."""
    return str(a).strip().upper() == b.upper()


def strip_unit(s):
    """Keep only digits, a decimal point and the letter 'E' (or 'e')."""
    return ''.join(ch for ch in s if ch.isdigit() or ch in {'.', 'E', 'e'})


def safeget(par, key, default=''):
    """
    Get a dictionary entry with a default.
    If the value is a string, strip single quotes.
    """
    val = par.get(key, default)
    if isinstance(val, str):
        val = val.replace("'", "")
    return val


# ----------------------------------------------------------------------
# Main routine â exactly the same logic as the MATLAB code
# ----------------------------------------------------------------------
def xeprpar(par):
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


# ----------------------------------------------------------------------
# Example usage â replace dummy_par with real data read from a .par file
# ----------------------------------------------------------------------
if __name__ == "__main__":
    dummy_par = {
        'MWFQ': '123.45',
        'CenterField': '3400 G',
        'IKKF': 'CPLX',
        'TITL': 'Test experiment',
        'RES': '512',
        'HCF': '3300',
        'HSW': '200',
        'JSS': '2',
    }

    out = xeprpar(dummy_par)
    for k, v in out.items():
        print(f"{k:15s}: {v}")



import os


def parse_dsc(file_path: str) -> dict:
    """
    Read the whole .DSC file into memory and parse it.

    Parameters
    ----------
    file_path : str
        Path to the .DSC file.

    Returns
    -------
    dict
        Dictionary of parameter names to their values.
    """
    # 1ï¸â£  Verify the file exists (os only)
    if not os.path.isfile(file_path):
        raise IOError(f"File not found: {file_path}")

    # 2ï¸â£  Load the entire file in one go
    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        raw_data = fh.read()

    # 3ï¸â£  Split into lines (Pythonâs builtâin str.splitlines)
    lines = raw_data.splitlines()

    # 4ï¸â£  Parse
    params = {}
    previous_key = None

    for raw_line in lines:
        line = raw_line.strip()
        if not line:                 # skip empty lines
            continue

        # Skip lines that start with a forbidden character
        if line[0] in _FORBIDDEN:
            continue
            
        # Split at the first whitespace
        parts = line.split(None, 1)
        if not parts:                # malformed line â ignore
            continue

        key, rest = parts[0], parts[1] if len(parts) > 1 else ""

        # New field â key does NOT start with a digit
        if not key[0].isdigit():
            previous_key = key
            val = rest.strip()
            val = _strip_quotes(val)
            params[key] = val
        # Continuation line â key starts with a digit
        else:
            if previous_key is None:  # no preceding key â ignore
                continue
            existing = params.get(previous_key, "")
            params[previous_key] = (existing + rest).strip()

    return params



# ----------------------------------------------------------------------
#  Demo â run the function from the command line
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python xeprdsc.py <file.DSC>")
        sys.exit(1)

    dsc_file = sys.argv[1]
    try:
        out = xeprdsc(dsc_file)
    except IOError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    # Prettyâprint the resulting dictionary (simple manual formatting)
    for k in sorted(out):
        print("{:<15} : {}".format(k, out[k]))

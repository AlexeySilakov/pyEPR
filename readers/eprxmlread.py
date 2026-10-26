"""Magnettech MS-5000 data XML (<ESRXmlFile>) reader, the companion of eprcsvread."""
import base64
import xml.etree.ElementTree as ET
import numpy as np
from scipy.signal import lfilter


def xmlread(fname):
    """[(ax, data, dsc), ...], one per <Measurement>, with ax/data/dsc shaped like eprcsvread.csvread.

    The file keeps time traces, not a field-swept spectrum: BField(t) and 'MWAbsorption (1st harm.)'(t),
    each with time = XOffset + XSlope*index and stored as little-endian doubles, base64-encoded
    one value per 12-character block. The absorption trace is mapped onto field through BField(t), then
    resampled (same number of points) onto a uniform field grid: the sweep rate is not constant, and
    the Hilbert transform used for phasing needs even spacing.
    Lock-in quadratures: the instrument's trace is S*cos(th0) + C*sin(th0) of the 'sinus'/'cosinus'
    traces (th0 = Measurement Phase), RC-filtered (FilterType DIGRC, tau = FilterPrm0), shifted back by tau
    and minus its starting value. ax['quad'] is the 90-degree partner C*cos(th0) - S*sin(th0), treated the
    same way and on the same grid, so data*cos(th - th0) + quad*sin(th - th0) is the spectrum at lock-in
    phase th; ax['phase0'] = th0 in degrees. quad is None when the file has no quadrature traces.
    dsc: 'Measurement.<attribute>' -> [value, ''] and 'Recipe.<parameter>' -> [value, unit]."""
    out = []
    for meas in ET.parse(fname).getroot().iter('Measurement'):
        dsc = {f"Measurement.{k}": [v, ''] for k, v in meas.attrib.items()}
        for p in meas.iter('Param'):
            dsc[f"Recipe.{p.get('Name')}"] = [(p.text or '').strip(), p.get('Unit') or '']
        curves = {}
        for c in meas.iter('Curve'):
            text = (c.text or '').strip()
            if not text or c.get('Compression', 'Base64') != 'Base64':
                continue
            v = np.frombuffer(b''.join(base64.b64decode(text[i:i + 12]) for i in range(0, len(text), 12)), '<f8')
            curves[c.get('YType')] = (float(c.get('XOffset', 0)) + float(c.get('XSlope', 0))*np.arange(v.size), v)
        if 'BField' not in curves or 'MW_Absorption' not in curves:
            raise ValueError(f"measurement '{meas.get('Name')}' has no BField / MW_Absorption trace")
        (tB, B), (tA, A) = curves['BField'], curves['MW_Absorption']
        phase0 = float(meas.get('Phase', 0) or 0)
        Q = None
        if ('MW_AbsorptionSinus' in curves and 'MW_AbsorptionCosinus' in curves
                and curves['MW_AbsorptionSinus'][1].size == curves['MW_AbsorptionCosinus'][1].size > 1):
            (tS, S), (_, C) = curves['MW_AbsorptionSinus'], curves['MW_AbsorptionCosinus']
            Q = C*np.cos(np.deg2rad(phase0)) - S*np.sin(np.deg2rad(phase0))
            tau = float(meas.get('FilterPrm0', 0) or 0) if meas.get('FilterType') == 'DIGRC' else 0.0
            if tau > 0:          # first-order RC, then the instrument's lag compensation
                a = (tS[1] - tS[0])/(tau + tS[1] - tS[0])
                Q = lfilter([a], [1, a - 1], Q, zi=[Q[0]*(1 - a)])[0]
                tS = tS - tau
            Q = np.interp(tA, tS, Q)
            Q -= Q[:20].mean()
        keep = (tA >= tB[0]) & (tA <= tB[-1])          # no extrapolation of the field trace
        BA, A = np.interp(tA[keep], tB, B), A[keep]
        order = np.argsort(BA, kind='stable')
        x = np.linspace(BA[order[0]], BA[order[-1]], A.size)
        ax = {'x': x, 'xlabel': f"BField [{meas.get('XUnit', 'mT')}]", 'ylabel': 'MW_Absorption []',
              'freq1': float(meas.get('MwFreq', 0) or 0), 'mwpower': 0.0, 'title': meas.get('Name', ''),
              'quad': None if Q is None else np.interp(x, BA[order], Q[keep][order]), 'phase0': phase0}
        try:
            ax['mwpower'] = float(dsc.get('Recipe.MicrowavePower', ['0'])[0])
        except ValueError:
            pass
        out.append((ax, np.interp(x, BA[order], A[order]), dsc))
    return out

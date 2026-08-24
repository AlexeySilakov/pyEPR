# -*- coding: utf-8 -*-
"""
Worker side of the grid search: everything one grid point needs, with no wx
anywhere in it.

Kept out of grid_search.py because multiprocessing on Windows uses 'spawn'
rather than fork
The single-process path runs these same functions in the GUI process, so the
two modes execute identical code and cannot drift apart.

@author: Alexey Silakov
"""
import copy
from multiprocessing import shared_memory

#import os
#import shutil
#import tempfile
#import time

import numpy as np

import mathfunctions as mf
import sys_functions as sysfun
from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim


class SilentSimError(Exception):
    pass

def silent_error_func(message):
    raise SilentSimError(message)

def log_gridpoint_error(error_log, where, exc):

    if error_log is None:
        return
    msg = f"{where}: {type(exc).__name__}: {exc}"
    error_log[msg] = error_log.get(msg, 0) + 1


# ---------------------------------------------------------------------------
# Sys-dict manipulation. These live here rather than on the panel so the
# worker can reach them; GridSearchFrame keeps its methods as one-line
# forwards, so every existing call site is unchanged.
# ---------------------------------------------------------------------------
def flatten_sys_params(sys_dict):
    """
    Flatten a Sys parameter dict (as stored in PropGridPanel.parameters,
    i.e. {'spin(1)': {...}, 'nuc(1)': {...}, ...}) into scannable scalar
    leaves. Only plain float/int values and elements of 1-D np.ndarray
    values are included; bools, strings and [value, choices] entries
    (e.g. 'Nucs', 'useFor') are skipped since a numeric grid scan over
    them isn't meaningful.
    """
    out = []
    for cat, sub in sys_dict.items():
        if not isinstance(sub, dict):
            continue
        for key, val in sub.items():
            if isinstance(val, bool):
                continue
            if isinstance(val, (int, float, np.floating, np.integer)):
                out.append({'kind': 'sysleaf', 'cat': cat, 'key': key, 'idx': None,
                            'value': float(val)})
            elif isinstance(val, np.ndarray) and val.ndim == 1:
                for ii in range(val.size):
                    out.append({'kind': 'sysleaf', 'cat': cat, 'key': key, 'idx': ii,
                                'value': float(val[ii])})
    return out


def param_leaf_label(p):
    return p['key'] if p['idx'] is None else f"{p['key']}({p['idx'] + 1})"


def qualified_label(p):
    return f"{p['cat']}.{param_leaf_label(p)}"


def apply_param_value(sys_dict, param, value):
    if param.get('kind', 'sysleaf') != 'sysleaf':
        raise ValueError(f"apply_param_value only applies to 'sysleaf' params, got {param}")
    sub = sys_dict[param['cat']]
    if param['idx'] is None:
        sub[param['key']] = float(value)
    else:
        arr = sub[param['key']]
        if not isinstance(arr, np.ndarray):
            arr = np.array(arr, dtype=float)
            sub[param['key']] = arr
        arr[param['idx']] = float(value)


def apply_functions_to_sys_dict(sys_dict, func_rows, variable_values, error_log=None):
    """
    If Function Modifiers are defined (func_rows non-empty), evaluate them
    -- using sys_functions.evaluate_functions(), the same UI-free engine
    classFunctionModPG.py uses -- and patch the results into `sys_dict` in
    place, resolving any "other property" references directly against
    `sys_dict` itself (so e.g. a formula referencing "spin(1).S" picks up
    whatever value a scanned sysleaf parameter already patched into
    sys_dict for this grid point).

    Returns True on success. Returns False (sys_dict left only partially
    patched) if any row failed to evaluate -- callers should treat that
    grid point as a failure (NaN) rather than silently simulating with
    stale/wrong values for the properties that were supposed to be driven.
    """
    if not func_rows:
        return True

    # qualified 'cat.leaf' label (e.g. 'nuc(1).A(1)') -> sysleaf descriptor,
    # for patching the results back in via apply_param_value below
    syspath_map = {qualified_label(p): p for p in flatten_sys_params(sys_dict)}
    known_labels = list(syspath_map.keys())

    def lookup(label):
        p = syspath_map.get(label)
        if p is None:
            raise sysfun.UnresolvedReference(label)
        sub = sys_dict[p['cat']]
        return float(sub[p['key']]) if p['idx'] is None else float(sub[p['key']][p['idx']])

    result = sysfun.evaluate_functions(func_rows, variable_values, lookup,
                                       known_labels=known_labels)

    if result['errors']:
        for idx, msg in result['errors'].items():
            lbl, expr = func_rows[idx]
            log_gridpoint_error(error_log, f'function "{lbl} = {expr}"', RuntimeError(msg))
        return False

    for label, value in result['label_values'].items():
        p = syspath_map.get(label)
        if p is not None:
            apply_param_value(sys_dict, p, value)
    return True


# ---------------------------------------------------------------------------
# Per-dataset records
# ---------------------------------------------------------------------------
def dataset_records(active_data, region_cache):
    """The slice of each dataset a grid point actually reads.

    Deliberately not the dataset dicts themselves: run_gridpoint only ever
    wants a handful of scalars, the orientation-selection grid and the
    pre-sliced region, while dd['data'] / dd['fftdata'] hold the full 2D maps
    the region was cut from. Those never need to cross a process boundary,
    and leaving them out is the difference between a payload of kilobytes and
    one of tens of megabytes per worker.
    """
    recs = []
    for dd, region in zip(active_data, region_cache):
        if dd['isfft']:
            npoints = len(dd['ax']['x'])
        elif dd.get('fftax'):
            npoints = len(dd['fftax']['x'])
        else:
            npoints = len(dd['ax']['x'])
        recs.append({
            'title': dd.get('title', dd.get('fname', '?')),
            'tau': dd['tau'],
            'field': dd['field'],
            'freq': dd['freq'],
            'fmax': dd['fmax'],
            'nPoints': npoints,
            'orisel': dd.get('orisel'),
            'region': region,
        })
    return recs


# ---------------------------------------------------------------------------
# Stored simulations
# ---------------------------------------------------------------------------
#def store_root():
#    """Scratch area for stored simulations, kept inside the project folder
#    rather than the system temp directory."""
#    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '_scan_cache')


#def write_point(dirpath, i, j, slot):
#    """Park one grid point's per-dataset simulations on disk. Returns the
#    number of bytes written, or 0 if the point produced nothing worth
#    keeping."""
#    present = np.array([s is not None for s in slot], dtype=bool)
#    if not present.any():
#        return 0
#    arrays = {'present': present}
#    for k, s in enumerate(slot):
#        if s is not None:
#            arrays[f'k{k}'] = np.asarray(s)
#    path = os.path.join(dirpath, f'p_{i}_{j}.npz')
#    np.savez(path, **arrays)
#    try:
#        return os.path.getsize(path)
#    except OSError:
#        return 0


class SimStore:
    """Every grid point's simulations in one shared block.

    The total is known before a scan starts -- n1*n2 points by the summed
    region size of every dataset -- so the block is allocated once and each
    (i, j, k) owns a fixed, disjoint slice. Workers write straight into their
    slice: nothing is pickled, nothing is copied through a pipe, nothing
    touches the disk, and no two points ever contend for the same bytes.

    One class, two roles. The panel constructs it with no `name` and owns the
    segment; a worker passes the name back in (see `layout`) and attaches to
    the same bytes. `owner` decides who is allowed to unlink it.

    Kept dict-shaped (len / bool / get / in) so the residual viewer cannot
    tell it from the plain {} it replaced.
    """

    # Display data only: these arrays feed the residual maps, while the fit
    # itself is already reduced to SumSqRes/n by the time one is stored.
    DTYPE = np.float32

    def __init__(self, shapes, n1, n2, name=None):
        self.shapes = [tuple(int(v) for v in s) for s in shapes]
        self.sizes = [int(np.prod(s)) for s in self.shapes]
        self.offsets = [0]
        for size in self.sizes:
            self.offsets.append(self.offsets[-1] + size)
        self.stride = self.offsets[-1]          # floats per grid point
        self.n1, self.n2 = int(n1), int(n2)

        count = self.n1 * self.n2 * self.stride
        self.nbytes = count * np.dtype(self.DTYPE).itemsize
        if name is None:
            # size=0 is rejected, and a scan with no region to store is not
            # worth a special case anywhere else.
            self.shm = shared_memory.SharedMemory(create=True,
                                                  size=max(self.nbytes, 1))
            self.owner = True
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.owner = False
        self.flat = np.ndarray(count, dtype=self.DTYPE, buffer=self.shm.buf)
        self._present = {}

    @property
    def layout(self):
        """Everything a worker needs to attach to this block, as a dict that
        feeds straight back into the constructor."""
        return {'shapes': self.shapes, 'n1': self.n1, 'n2': self.n2,
                'name': self.shm.name}

    def view(self, i, j, k):
        """Dataset k's slice at grid point (i, j), as a live array."""
        base = (i * self.n2 + j) * self.stride + self.offsets[k]
        return self.flat[base:base + self.sizes[k]].reshape(self.shapes[k])

    def note(self, i, j, present):
        """Record which datasets a worker actually wrote at this point. The
        block itself carries no presence flag, and reading a slice nobody
        wrote would hand back zeros as if they were a simulation."""
        if any(present):
            self._present[(i, j)] = tuple(present)

    def get(self, key, default=None):
        """Copies, not views: the residual code is free to work in place, and
        an outstanding view would block close() when the store is disposed.
        One grid point is a few hundred kB, so the copy is not worth avoiding.
        """
        present = self._present.get(key)
        if present is None:
            return default
        i, j = key
        return [self.view(i, j, k).copy() if p else None
                for k, p in enumerate(present)]

    def keys(self):
        return self._present.keys()

    def __contains__(self, key):
        return key in self._present

    def __len__(self):
        return len(self._present)

    def __bool__(self):
        return bool(self._present)

    def dispose(self):
        """Release the block. The view has to go first -- close() refuses
        while any exported pointer into the buffer is still alive."""
        self._present = {}
        self.flat = None
        try:
            self.shm.close()
            if self.owner:
                self.shm.unlink()   # no-op on Windows, required on POSIX
        except (BufferError, FileNotFoundError, OSError):
            pass
#def sweep_store_root(keep=None, stale_age=86400.0, include_own=True):
#    """Remove leftover run directories. A hard cancel terminates workers
#    mid-write and a crash leaves the whole directory behind, so something has
#    to clear them out, and scan start is the natural moment.
#
#    Anything more than a day old is debris from a dead session and always
#    goes. Our own earlier runs go too unless include_own is False, which is
#    what a closing panel wants: it disposes the stores it created by name, and
#    a second Grid Search window in the same process has the same pid but
#    stores of its own that must survive.
#
#    A recent directory belonging to another pid is always left alone -- it may
#    be a second copy of the program running right now, and deleting its store
#    out from under it would break its residual maps.
#    """
#    base = store_root()
#    if not os.path.isdir(base):
#        return
#    mine = f'run_{os.getpid()}_'
#    now = time.time()
#    for name in os.listdir(base):
#        path = os.path.join(base, name)
#        if path == keep or not name.startswith('run_'):
#            continue
#        try:
#            stale = (now - os.path.getmtime(path)) > stale_age
#        except OSError:
#            stale = False
#        if (include_own and name.startswith(mine)) or stale:
#            shutil.rmtree(path, ignore_errors=True)


# ---------------------------------------------------------------------------
# One grid point
# ---------------------------------------------------------------------------
def run_gridpoint(hs, Sys, exp_template, opt_template, records, error_log=None,
                  sim_store=None, stats_out=None):
    total_SumSqRes = 0.0
    total_n = 0
    norm_acc = 0.0    # sum of N_k * (SumSqRes_k / null_k)
    hs.Sys = Sys
    try:
        hs.preCompute()
    except Exception as e:
        log_gridpoint_error(error_log, "preCompute", e)
        return np.nan, 0, np.nan

    for k, rec in enumerate(records):
        region = rec['region']
        if region is None:
            continue
        Exp = expPar()
        Exp.setDict(copy.deepcopy(exp_template))
        if Exp.tau < 0:
            Exp.tau = rec['tau']
        if Exp.Field < 0:
            Exp.Field = rec['field']
        if Exp.mwFreq < 0:
            Exp.mwFreq = rec['freq']
        if Exp.MaxFreq < 0:
            Exp.MaxFreq = rec['fmax']
        if Exp.nPoints < 0:
            Exp.nPoints = rec['nPoints']

        Opt = optHYSCORE()
        Opt.setFromCtrl(copy.deepcopy(opt_template))
        if rec.get('orisel') is not None:
            Opt.OriSelInp = rec['orisel']

        hs.Exp = Exp
        hs.Opt = Opt

        ### Do run
        try:
            hs.reRun()
        except Exception as e:
            log_gridpoint_error(error_log, f"reRun ({rec['title']})", e)
            continue

        if hs.Spectrum is None:
            continue

        # Calculate RMSD
        mask = region['mask']
        exp_region = region['exp']
        offset = region['offset']

        sim_on_exp = mf.adaptSim(region['x'], region['y'], hs.Spectrum, hs.X, hs.Y)
        SumSqRes, n, scaled_sim = mf.calcSSR(exp_region, sim_on_exp, mask, offset)

        if n > 0:
            total_SumSqRes += SumSqRes
            total_n += n
            if sim_store is not None:
                sim_store[k] = scaled_sim
            if stats_out is not None:
                stats_out[k] = (SumSqRes, n)
            norm_acc += n * (SumSqRes / region['null'])

    if total_n == 0:
        return np.nan, 0, np.nan
    rmsd_norm = float(np.sqrt(norm_acc / total_n))
    rmsd = float(np.sqrt(total_SumSqRes / total_n))
    # total_n - number of points accross all datasets
    return rmsd, total_n, rmsd_norm


# ---------------------------------------------------------------------------
# Pool plumbing
# ---------------------------------------------------------------------------
# Per-process state. In a spawned worker this is filled once by init_worker
# and reused for every point that worker is handed, so the HYSCOREsim
# instance -- and with it the grid cache ensure_grid() maintains -- survives
# between points exactly as it did in the old single-threaded loop.
_CTX = {}


def init_worker(payload):
    hs = HYSCOREsim(errorFunc=silent_error_func)
    hs.verbose = False
    _CTX.clear()
    _CTX.update(payload)
    _CTX['hs'] = hs
    layout = payload.get('store_layout')
    # The SimStore has to be held, not just used: dropping it closes the
    # segment and leaves the numpy view pointing at freed memory.
    _CTX['block'] = SimStore(**layout) if layout else None

def release_worker():
    """Let go of the shared block. Pool workers exit and take their
    attachment with them; the single-process path runs inside the GUI, which
    does not, so it has to release by hand or the segment outlives the scan."""
    block = _CTX.pop('block', None)
    if block is not None:
        block.dispose()
    _CTX.clear()
    
def run_point(task):
    """One grid point, start to finish.

    Returns a tuple rather than writing into shared arrays: the coordinator
    owns the result arrays and is the only thing that touches them, so there
    is nothing to synchronise. Errors come back as a small {message: count}
    dict for the coordinator to merge and print.
    """
    i, j, v1, v2 = task
    hs = _CTX['hs']
    records = _CTX['records']
    store_dir = _CTX.get('store_dir')
    error_log = {}

    sys_dict = copy.deepcopy(_CTX['sys_template'])
    variable_values = dict(_CTX['base_variable_values'])
    for p, v in ((_CTX['param1'], v1), (_CTX['param2'], v2)):
        if p.get('kind', 'sysleaf') == 'funcvar':
            variable_values[p['key']] = v
        else:
            apply_param_value(sys_dict, p, v)

    slot = [None] * len(records) if _CTX.get('block') is not None else None
    stats = [None] * len(records)

    func_ok = False
    try:
        func_ok = apply_functions_to_sys_dict(
            sys_dict, _CTX['func_rows'], variable_values, error_log=error_log)
        Sys = sysPar()
        Sys.setFromCtrl(sys_dict)
    except Exception as e:
        log_gridpoint_error(error_log, "run_point->Sys.setFromCtrl", e)

    t_rmsd, t_total_n, t_rmsd_norm = np.nan, 0, np.nan
    if func_ok:
        try:
            t_rmsd, t_total_n, t_rmsd_norm = run_gridpoint(
                hs, Sys, _CTX['exp_template'], _CTX['opt_template'], records,
                error_log=error_log, sim_store=slot, stats_out=stats)
        except Exception as e:
            log_gridpoint_error(error_log, "run_point->run_gridpoint", e)


    present = [False] * len(records)
    if slot is not None:
        block = _CTX['block']
        for k, sim in enumerate(slot):
            if sim is None:
                continue
            try:
                # float64 -> float32 happens on assignment into the view
                block.view(i, j, k)[...] = sim
                present[k] = True
            except Exception as e:
                log_gridpoint_error(error_log, "run_point->store", e)

    return (i, j, t_rmsd, t_total_n, t_rmsd_norm, stats, present, error_log)

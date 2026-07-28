"""
Pure (no wx dependency) expression-evaluation engine shared by the Sys
"Function Modifiers" panel (classFunctionModPG.py, live wx property grid)
and the Grid Search worker thread (grid_search.py, plain dicts). Kept
UI-free so it's safe to call from a background thread.

An expression drives a target property's value from a mix of custom scalar
parameters (e.g. r, a) and/or other property values (e.g. "spin(1).S"),
referenced by their qualified "<category>.<leaf>" label. Multiple driven
targets are evaluated in dependency order, with duplicate-target and
circular-reference detection.
"""
import math
import re

SAFE_FUNCS = {
    'sqrt': math.sqrt, 'exp': math.exp, 'log': math.log, 'log10': math.log10,
    'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
    'asin': math.asin, 'acos': math.acos, 'atan': math.atan, 'atan2': math.atan2,
    'pi': math.pi, 'e': math.e, 'abs': abs, 'min': min, 'max': max, 'pow': pow,
}


class UnresolvedReference(LookupError):
    """Raise this from a prop_value_lookup callable passed to
    evaluate_functions() when a referenced label can't be resolved."""
    pass


def extract_refs_and_freevars(expr, known_labels, safe_funcs=None):
    """Split `expr` into a form eval() can use plus the set of "free"
    (custom-parameter) names it references. Any literal occurrence of a
    string in `known_labels` (property labels, which may contain '(', ')'
    or '.' and so can't be picked up by a plain identifier regex) is
    replaced with an opaque placeholder token first, longest label first
    so one label can't accidentally match inside a longer one."""
    safe_funcs = safe_funcs if safe_funcs is not None else SAFE_FUNCS
    work = expr.split('=')[-1].strip().replace('^', '**')
    token_map = {}   # placeholder token -> qualified label
    for i, lbl in enumerate(sorted(set(known_labels), key=len, reverse=True)):
        if lbl and lbl in work:
            token = f"__PGREF{i}__"
            work = work.replace(lbl, token)
            token_map[token] = lbl
    freevars = set(re.findall(r"[A-Za-z_]\w*", work))
    freevars -= set(token_map.keys())
    freevars -= set(safe_funcs.keys())
    return work, token_map, freevars


def topo_sort(deps):
    """Kahn's algorithm. deps: node -> set of nodes that must be evaluated
    first. Returns (order, cyclic_nodes) -- cyclic_nodes lists every node
    that could not be scheduled because it's part of a cycle."""
    nodes = list(deps.keys())
    indeg = {n: len(deps[n]) for n in nodes}
    children = {n: [] for n in nodes}
    for n, ds in deps.items():
        for d in ds:
            if d in children:
                children[d].append(n)
    queue = [n for n in nodes if indeg[n] == 0]
    order = []
    while queue:
        n = queue.pop(0)
        order.append(n)
        for c in children[n]:
            indeg[c] -= 1
            if indeg[c] == 0:
                queue.append(c)
    cyclic = [n for n in nodes if n not in order]
    return order, cyclic


def evaluate_functions(rows, variable_values, prop_value_lookup,
                        known_labels=None, safe_funcs=None):
    """
    rows: ordered list of (target_label, expression_string) tuples (one
        per function-modifier row, in UI order -- may contain duplicate
        target_labels, which are detected and rejected).
    variable_values: dict of custom-parameter name -> float.
    prop_value_lookup: callable(label) -> float, used to resolve a
        reference to a property that ISN'T itself one of `rows`'s targets
        (i.e. a "base" value). Must raise LookupError (or the
        UnresolvedReference subclass) if `label` can't be resolved.
    known_labels: iterable of every property label that may legally be
        referenced from an expression. Defaults to just the targets in
        `rows`, which is enough for expressions that only reference other
        driven targets, but you'll usually want to also pass every plain
        property label so cross-references to non-driven properties can be
        recognized by extract_refs_and_freevars().

    Returns a dict:
      'values':       {row_index: float}   -- rows that evaluated OK
      'errors':       {row_index: str}     -- rows that failed, and why
      'label_values': {target_label: float} -- convenience map of OK rows
      'order':        [row_index, ...] evaluation order actually used
    """
    safe_funcs = safe_funcs if safe_funcs is not None else SAFE_FUNCS
    all_known = set(known_labels) if known_labels is not None else set()
    all_known |= {lbl for lbl, _expr in rows}

    label_counts = {}
    for lbl, _expr in rows:
        label_counts[lbl] = label_counts.get(lbl, 0) + 1
    duplicate_labels = {lbl for lbl, n in label_counts.items() if n > 1}

    parsed = {}   # row_idx -> (label, work, token_map, freevars)
    errors = {}
    for idx, (lbl, expr) in enumerate(rows):
        if lbl in duplicate_labels:
            errors[idx] = f'"{lbl}" is targeted by more than one function'
            continue
        work, token_map, freevars = extract_refs_and_freevars(expr, all_known, safe_funcs)
        parsed[idx] = (lbl, work, token_map, freevars)

    label_to_idx = {lbl: idx for idx, (lbl, *_r) in parsed.items()}
    deps = {idx: set() for idx in parsed}
    for idx, (lbl, _work, token_map, _freevars) in parsed.items():
        for _tok, ref_lbl in token_map.items():
            if ref_lbl in label_to_idx and ref_lbl != lbl:
                deps[idx].add(label_to_idx[ref_lbl])

    order, cyclic = topo_sort(deps)
    for idx in cyclic:
        errors[idx] = "circular reference between functions"

    values = {}
    label_values = {}
    for idx in order:
        lbl, work, token_map, freevars = parsed[idx]

        env = dict(safe_funcs)
        missing_var = None
        for name in freevars:
            if name in variable_values:
                env[name] = variable_values[name]
            else:
                missing_var = name
                break
        if missing_var is not None:
            errors[idx] = f'undefined parameter "{missing_var}"'
            continue

        unresolved = None
        for tok, ref_lbl in token_map.items():
            if ref_lbl in label_values:
                env[tok] = label_values[ref_lbl]
            else:
                try:
                    env[tok] = prop_value_lookup(ref_lbl)
                except LookupError:
                    unresolved = ref_lbl
                    break
        if unresolved is not None:
            errors[idx] = f'reference "{unresolved}" could not be resolved'
            continue

        try:
            result = float(eval(work, {"__builtins__": {}}, env))
        except Exception as e:
            errors[idx] = f"{type(e).__name__}: {e}"
            continue

        values[idx] = result
        label_values[lbl] = result

    return {'values': values, 'errors': errors, 'label_values': label_values, 'order': order}

# pyCWEPR

A wxPython program for loading, displaying and comparing CW EPR spectra, with
background subtraction, manual absorption/dispersion phase correction and
double integration. It uses the same layout and widgets as pyHYSCORE.

Start it with:

```
python pyCWEPR.py
```

Requirements: `pip install wxPython numpy scipy matplotlib`

## Folder contents

| File | Contents |
|---|---|
| `pyCWEPR.py` | The program: `MainFrame`, `FileBrowserPanel`, `MatplotlibPanel`, `TabulatedPanel`, processing (`MainFrame.process`), session XML |
| `theme/` | Shared with pyHYSCORE: `theme.py` (pill widgets, `Card`, `draw_curve_icon`, spin-box Enter handling), `classPropGridPanel.py`, `classFunctionModPG.py` |
| `spin_system/sys_functions.py` | Expression engine used by `classFunctionModPG` |
| `readers/` | `brukerread.py` (Bruker BES3T/ESP), `eprcsvread.py` (Magnettech CSV), `eprxmlread.py` (Magnettech data XML: `xmlread`, same `ax, data, dsc` shape as `csvread`) |
| `icons/` | File-list icons, found relative to the script (any working directory works) |

See the module map in [`README.md`](README.md) for the full folder layout.

## Window layout

```
┌ toolbar: [folder][up] | [save] | [deriv][abs][DI] …… [B/g][grid][theme] ┐
│◂│ file browser │ plot (checkboxes, canvas, nav toolbar) │ tabs │▸│
└ status: [message] [cursor B / g / y readout] ┘
```

## Input formats

* **Magnettech MS-5000 `*_result.csv`**: the field axis is in mT and the MW
  frequency comes from `Additional.Frequency`.
* **Plain CSV / TXT** (any `.csv` that is not a Magnettech result, or `.txt`,
  `.dat`, `.asc`). Columns are separated by `,`, `;`, tab or spaces; with `;`
  or tab a decimal comma also works.
  * `#` lines are metadata, read as `# key: value` and shown by **Header**. An
    optional row of column names comes next, then the numbers.
  * The first column is the field, in mT, or in G if its name says so. The
    spectrum is the column whose name contains "deriv", or else the second
    column. `# MW frequency [GHz]: …` sets the frequency.
  * This is the format **Export ASCII** writes, so exported files load back.
* **Bruker `.DSC/.DTA`** (and ESP `.par/.spc` through `brukerread`): a field
  axis in G is converted to mT. Complex data uses the real channel.
* **XML**
  * A *Magnettech data XML* (`<ESRXmlFile>`, saved next to the `_result.csv`)
    holds time traces, not a spectrum. The field trace B(t) maps the
    "MWAbsorption (1st harm.)" trace onto field. That is resampled, with the
    same number of points, onto a uniform field grid, because the sweep rate is
    not constant and the Hilbert phase step needs even spacing. The result is
    the CSV spectrum at about 3× its resolution and over the slightly wider
    measured field range, without the CSV's constant offset. The header popup
    shows the measurement attributes and recipe parameters with units. A file
    with several measurements adds one spectrum each.
  * The file also keeps the lock-in's raw sin/cos traces S and C. The
    instrument's trace is S·cos θ₀ + C·sin θ₀ (θ₀ = the file's `Phase`),
    RC-filtered (τ = `FilterPrm0`), shifted back by τ, minus its starting value.
    The reader builds the 90° partner Q = C·cos θ₀ − S·sin θ₀ the same way. Then
    `data·cos(θ−θ₀) + Q·sin(θ−θ₀)` is the spectrum at any lock-in phase θ, and
    θ = θ₀ gives exactly the instrument's trace.
  * A *pyCWEPR session* restores everything: spectra, backgrounds, processing
    settings and view.
  * Any other session in the same encoding (pyEPR, pyHYSCORE) adds its 1D
    spectra (its `Data` list) as new spectra, or as backgrounds with **+ BG**.

Ctrl- or Shift-click selects several files in the browser. **Load** and **+ BG**
then load all of them.

**Header** (in the file browser) shows a file's parameters in a popup without
loading the file. **Header** on the Data and Background tabs does the same for
loaded spectra.

## Processing chain (per spectrum, `MainFrame.process`)

0. **Lock-in phase** (Magnettech XML spectra only): `y = data·cos(θ−θ₀) + Q·sin(θ−θ₀)`.
1. **Background subtraction.** Each entry in the spectrum's background list has
   a field axis that is mapped as `B·ratio + shift`. The mapped curve is
   interpolated onto the spectrum's own field axis with a cubic spline, then
   `scale × background` is subtracted.
   * `ratio = −1` means use `f(spectrum)/f(background)`, so the same g-value
     lands on the same point. Any positive value overrides it.
   * Outside the background's field range the background counts as 0, and the
     status bar shows a warning.
2. **Derivative baseline.** A polynomial of the chosen order is fitted to the
   points inside the "Der.BL" regions and subtracted. Order −1 turns this off.
3. **Phase.** `y = Re[e^{iφ}(s + i·H[s])] = s·cos φ − H[s]·sin φ`, where H is the
   Hilbert transform of the corrected derivative `s`. `s` is zero-padded to 2N
   so the two ends of the spectrum do not wrap into each other.
4. **First integral** over the integration range gives the absorption.
5. **Absorption baseline.** A polynomial is fitted in the "Abs.BL" regions and
   subtracted.
6. **Second integral** gives the double-integral curve.
7. **Double-integral baseline.** A polynomial is fitted to the double integral
   in its own regions and subtracted. **DI** is the last value of the result,
   in data units × mT².

Defaults for a new spectrum:
* Derivative baseline: order 0, regions at the outer 10 % of each end of the
  field range.
* Absorption baseline: order 1, the same two regions.
* Double-integral baseline: order 0, the left 10 % region only.
* Integration range: the whole sweep.

Sessions saved before the double-integral baseline existed get its defaults
when loaded, which can shift their DI values slightly.

## Tabs

* **Spectrum** (above the tabs): selects the *current* spectrum. The
  Background, Phase and Integrate tabs edit this spectrum. Clicking a spectrum
  in the Data grid also selects it.
* In the Phase and Integrate tabs each group of controls sits on its own
  shaded card. **Right-click** a phase or integration-range number box to
  choose its arrow step.
* **Data**: Show, MW frequency, display Scale, display Y shift and a
  read-only DI for each spectrum. Also: delete, header, and reorder (▲▼).
  * **Compare** sets the normalisation: As is / Max (derivative peak-to-peak
    = 1) / DI (double integral = 1).
  * **Stack offset**: vertical step between overlaid spectra.
  * **Copy DI table** puts a tab-separated table on the clipboard.
  * **Export ASCII** saves the current spectrum's traces to a comma-separated
    `.csv`. You choose any of: raw, background, derivative (fully processed),
    absorption (after its baseline) and double integral.
    * `#` lines at the top record the title, source, frequency, phases,
      backgrounds, baseline settings, integration range and DI.
    * Next comes one row of column names, then `B [mT]` and the chosen traces.
    * Absorption and double integral exist only inside the integration range,
      so their cells are empty outside it.
    * Opens in Excel and Origin; with pandas use `read_csv(f, comment='#')`,
      with numpy `genfromtxt(f, delimiter=',', names=True)`.
    * Loading it back takes the derivative column as the spectrum.
* **Background**: the loaded background pool, and the list of backgrounds
  subtracted from the current spectrum (any number). Each entry has its own
  scale, frequency ratio and field shift. **Apply to all spectra** copies the
  list to every spectrum.
* **Phase**: two blocks, each with a slider (live update, 0.1° steps) and a
  spin box.
  * **Lock-in (modulation) phase θ**: absolute, starting at the file's value
    (shown next to it, with ↺ to go back). It is enabled only for spectra that
    carry sin/cos traces.
  * **Absorption / dispersion phase φ**: the Hilbert mixing. Watch the
    Absorption row while adjusting.
  * **Apply to all spectra** copies both. θ is ignored for spectra without
    traces.
* **Integrate**: the integration range, and polynomial order plus regions for
  the derivative, absorption and double-integral baselines.
  * Each fit region has its own **from / to** row in mT, with × to delete it.
    **+ Add range** adds one after the last region, 5 % of the field span
    wide. **Clear all** removes them all. Regions can also be edited on the
    plot (see below), and the rows follow.
  * The tab scrolls when there are many regions.
  * The DI readout is shown here. **Apply to all spectra** copies the
    integration settings to every spectrum.

* **Math**: an expression over loaded spectra makes a new, *live* spectrum.
  The syntax help is the expression box's tooltip.
  * `s1`, `s2`, … are spectra by their number in the Data list. Any other name
    (`a`, `k`, `c1`, …) is an adjustable constant.
  * Also allowed: numbers, `+ - * / **`, `abs max min mean sqrt`, and
    `diff`, the derivative d/dB along the result's field axis (per mT). The
    expression is checked before it is compiled, so nothing else can run.
    Examples: `s1 - a*s2`, `(s1 + s2)/2`, `k*(s3 - s4)`, `diff(s1)`.
  * The pill buttons under the box come in four rows of four: `s1`–`s4`,
    `a`–`d`, `abs() diff() ( )`, and `+ - / *`. Each inserts its text at the
    cursor. Functions leave the cursor between the brackets.
  * **Average** and **Weighted** fill in an expression over every loaded
    spectrum, whether shown or not, excluding math spectra. Weighted gives
    `c1*s1 + c2*s2 + …`, with one slider per weight.
  * **Create** adds the expression as a new math spectrum. It takes the field
    axis and frequency of the first spectrum in the expression.
  * **Update** (or Enter) applies the expression to the current spectrum if
    that is a math spectrum, otherwise to the last math spectrum in the Data
    list. If no math spectrum exists yet, it creates one.
  * Every other source is cubic-spline interpolated onto that axis, and is 0
    outside its own range, with a warning.
    * **Align Field** matches points in mT.
    * **Align g** first scales each source's field by f(result)/f(source).
  * **Stage Raw** uses the loaded data, after the lock-in phase. **Stage
    Processed** uses each source's final derivative, after its own background,
    baseline and phase.
  * While a math spectrum is current, the tab shows its sources (`s1 = title`)
    and one row per constant. Each row has a slider, a value box (right-click
    it to set the arrow step, as on the other spin boxes) and a "± range" that
    sets the slider's span. The plot follows as you drag.
    * Editing the expression and pressing Enter changes it.
    * Align and Stage switch it live.
    * It recomputes on every redraw, so changing a source's processing
      (Stage Processed) updates it too.
  * Sources are remembered by title, so reordering the Data list doesn't
    break a math spectrum. The expression is shown renumbered.
  * If a source is deleted, the last result stays and a warning appears.
  * A math spectrum gets its own background, phase and integration settings,
    and is saved in sessions. Math spectra cannot be sources of other math
    spectra.

## Plot

* Rows: Derivative, Absorption and Double integral. Each is switched by its
  own toggle on the main toolbar. Each icon shows the row's line shape: a
  derivative wiggle, an absorption peak, an integral step. Any combination
  works, but at least one row stays on.
* **X axis**: the **B / g** toggle on the toolbar switches between field (B,
  off) and g (on). The g axis uses each spectrum's own MW frequency, so
  spectra recorded at different frequencies line up.
* All shown spectra are overlaid. The current spectrum is drawn thicker.
* **Raw & BG** adds the raw spectrum and the summed background.
* **Baselines** adds the curves before baseline subtraction, the fitted
  polynomials (dotted) and the fit regions (grey).
* Dashed vertical lines mark the integration range.
* In the Double integral row, a dash-dot line marks the level at which the
  current spectrum's DI is read (the last point of the corrected curve), with
  "DI = …" written at its right end. The label gives the raw DI, whatever the
  Compare normalisation.
* Editing on the plot (current spectrum, Zoom/Pan off; the cursor changes over
  anything you can grab):
  * Drag a grey baseline region to move it, or drag its edge to resize it.
    Regions are shown, and editable, while **Baselines** is checked:
    derivative regions in the top row, absorption regions in the Absorption
    row, double-integral regions in the Double integral row.
  * **Shift-drag** on empty space in one of those rows adds a region there.
  * **Right-click** a region to delete it. If a spectrum line is under the
    cursor, a menu offers the line's colour and "Delete baseline region".
  * Drag a dashed line to move that integration limit. A limit wins over a
    region edge at the same place.
* **Double-click a spectrum line**, in any row, to make it the current
  spectrum. Where lines overlap, a non-current one is chosen. **Double-click
  empty space** to go back to the Home view, like the plot toolbar's Home.
* **Ctrl-Z** switches on zoom-to-rectangle for one zoom. After that zoom the
  mouse goes back to normal. In a text box, Ctrl-Z stays undo.
* **Drag a spectrum line up or down** (any spectrum, any row; the cursor shows
  ↕) to set its **Y shift**, a display-only vertical offset that is also in
  the Data grid.
  * It is measured in units of the tallest trace of each row, so one value
    moves the spectrum by the same proportion in all three rows.
  * It adds to the Stack offset. A line beats the inside of a baseline region,
    but region edges and integration limits still come first.
  * Hold Shift to draw a new region across a line instead.
* **Right-click a spectrum line** to pick its colour, or go back to the default.
  The colour is saved in the session.
* An x-zoom is kept while you adjust parameters, and y re-scales to the
  visible window. Home returns to the full view.
* **Right-click an axis** (its tick labels or title) to open a small dialog
  with Min / Max and **Autoscale**.
  * The x axis is shared: its range applies to all rows. It is kept
    separately for the field and g axes.
  * A y range applies to that row only. Rows without a fixed y range keep
    autoscaling to the visible x window.
  * A fixed range survives redraws and is saved in sessions. Tick Autoscale
    to release it.

## Not (yet) included

* Automatic phase or background-scale fitting. All corrections are manual.
* DI normalisation by acquisition parameters (modulation, power, Q, gain,
  temperature). Only the raw DI is reported, as requested.

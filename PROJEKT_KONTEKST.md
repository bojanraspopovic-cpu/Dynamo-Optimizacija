# DTM Smoothing — Projektni Kontekst

## Ko radi i kako

Građevinski inženjer infrastrukture (3-6 god. iskustva), Strabag Beč.
Alati: Civil 3D, AutoCAD, ArcGIS, QGIS, Dynamo, Python osnove.

**Kada pomažeš s kodom:** Ne objašnjavaj od nule. Idi direktno na rješenje, objasni samo ključne dijelove. Poznaje C3D, koridore i infrastrukturno projektovanje.

---

## Projektni zadatak

Klijent geodetski snimio trasu autoputa 6km. Snimak: svakih 5m duž trase, 10-15 tačaka po poprečnom presjeku, na mjestima saobraćajnih traka. Rezultat: as-built DTM — istrošen asfalt sa plastičnom deformacijom (bez velikih oštećenja).

Cilj: predložiti novi glatki DTM koji definišu geometriju overlay sloja.

---

## Geometrija trake

Radimo po traci odvojeno. Kontekst: **desna traka autoputa**.

- **Pivot (unutrašnja ivica)** — lijeva ivica desne trake (sredina kolnika). Rotacijska tačka od koje se baca poprečni nagib.
- **Spoljna ivica** — desna ivica desne trake.
- Poprečni nagib je uvijek upravan na tangentu trase u toj tački (ne u globalnom X/Y smjeru).
- Širina trake nije konstantna — varira na krivinama i prijelaznicama.

**Signed nagib:**
```
nagib = (Z_spoljna - Z_pivot) / sirina
Pozitivan: Z pada od pivota prema van (konveksna krivina)
Negativan: Z pada od van prema pivotu (konkavna krivina)
```

**Odvodnja check:** `abs(nagib) >= 0.025` (2.5%)

Na PRAVAC/KRIVINA dionicama: ako je `abs(nagib_smooth) < 2.5%` → automatski postavi nagib na `sign(nagib) × 2.5%`, prilagodi spoljnu Z, provjeri constraint → status `AUTO-ADJUSTED`.

---

## Constraint

`[-1cm / +3cm]` u odnosu na originalne Z vrijednosti.

Provjerava se na:
1. Pivot Z nakon smoothinga
2. Spoljna Z

Tačke između ivica ne trebaju zasebnu provjeru — automatski prolaze ako ivice prolaze.

---

## Ulazni podaci

### Skripta 0 (pivot smoothing):
```
pivot_raw.csv      → originalna 3D pivot polilinija iz C3D (X, Y, Z)
pivot_minus1.csv   → pivot linija offsetovana -1cm (donja granica koridora)
pivot_plus3.csv    → pivot linija offsetovana +3cm (gornja granica koridora)
```
> `vertikalni_elementi.csv` se NE koristi. Skripta sama detektuje geometriju iz podataka.

### Skripta 1 (glavna DTM skripta):
```
tacke.csv          → sve XYZ tačke trake (geodetski snimak)
pivot_smooth.csv   → output Skripte 0
spoljna.csv        → 3D polilinija spoljne ivice iz C3D (X, Y, Z)
alignment_info.csv → od_stat, do_stat, tip, smjer_nagiba
                     tip:          PRAVAC | KRIVINA | PRIJELAZNICA
                     smjer_nagiba:  1 = od pivota prema van (konveksna)
                                   -1 = od van prema pivotu (konkavna)
                                    0 = prijelaznica
```

**Format svih fajlova:** X,Y,Z po redu (sa ili bez zaglavlja).

**Stacionaže se mjere duž pivot 3D polilinije.**

---

## Output

### Skripta 0:
```
pivot_smooth.csv  → optimizovana pivot linija (X, Y, novi Z)
pivot_qa.csv      → po tački: orig_Z | novi_Z | delta |
                    status: NORMALNO | OUTLIER | KORIDOR_GRANICA
```

### Skripta 1:
```
novi_dtm.csv      → originalne XY, nove Z za sve tačke
qa_report.csv     → po presjeku:
                    stacionaža | tip_dionice | nagib_stari | nagib_novi |
                    delta_pivot | delta_spoljna | odvodnja_ok | status
                    status:      OK | AUTO-ADJUSTED | MANUAL CHECK | PRIJELAZNICA
                    odvodnja_ok: DA | NE | N/A
```

---

## SKRIPTA 0 — Pivot Smoothing

### Problem
Originalna pivot linija ima blagi šum (normalno) i mjestimično ekstremne anomalije — "zubi ajkule" koji nisu stvarna geometrija (vegetacija, loš snimak).

### Tri situacije:
- **Situacija A — Normalna:** blagi šum oko trenda → smooth prati trend
- **Situacija B — Outlier:** nagla skokovi koji nisu geometrija → ignoriši Z, provuci kroz koridor
- **Situacija C — Vertikalni luk:** Z prati krivinu → smooth mora pratiti krivinu, ne ispravljati je u pravac

### Outlier detekcija (bez external alignment info):
Koristiti **rolling median** na prozoru ±10 tačaka:
```
lokalni_median = rolling_median(Z, window=21)
odstupanje = abs(Z_i - lokalni_median)
ako odstupanje > prag (npr. 5cm) → OUTLIER
```
Rolling median je robustan — vertikalni luk ne uzrokuje lažne flagove jer median prati krivinu.

Za outlier tačke: koridor se bazira na **interpolaciji iz susjednih normalnih tačaka**, ne na izmjerenoj Z.

### Matematički pristup — QP (cvxpy + OSQP):
```
Minimizuj:  Σ (Δ²Z_i)²   ← minimizacija krivine (druga diferencija)
Uz uslov:   lb_i ≤ Z_i ≤ ub_i   za svaku tačku i

Normalna tačka:  lb = Z_orig - 1cm,  ub = Z_orig + 3cm
Outlier tačka:   lb = Z_interpol - 1cm,  ub = Z_interpol + 3cm
```
`cvxpy + OSQP` solver — hard constraints su nativni, nema rizika prekoračenja koridora.

> **Ne koristiti:** UnivariateSpline (ne podržava hard bounds nativno).
> **Ne koristiti:** KDTree (nije potreban za ovu veličinu podataka, ~6000 tačaka).

---

## SKRIPTA 1 — Glavna DTM Skripta

### Processing pipeline:

**1.** Učitaj tacke.csv, pivot_smooth.csv, spoljna.csv, alignment_info.csv

**2.** Reparametrizuj obje polilinije na jednake intervale (1m):
```python
def reparametrizuj(polilinija, interval=1.0):
    duzine = np.cumsum(np.r_[0, np.linalg.norm(
        np.diff(polilinija[:,:2], axis=0), axis=1)])
    ukupna = duzine[-1]
    nove_duzine = np.arange(0, ukupna, interval)
    x = interp1d(duzine, polilinija[:,0])(nove_duzine)
    y = interp1d(duzine, polilinija[:,1])(nove_duzine)
    z = interp1d(duzine, polilinija[:,2])(nove_duzine)
    return np.column_stack([x, y, z])
```

**3.** Provjeri orijentaciju normale jednom na početku — spoljna ivica mora biti konzistentno s jedne strane.

**4.** Za svaku stacionažu (na svakih 5m duž pivot polilinije):
- a. Tangenta — centrirana diferencija sa širim prozorom `(p[i+3] - p[i-3])`
- b. Normala na tangentu
- c. Stvarna širina = udaljenost do spoljne polilinije u smjeru normale
- d. Tip dionice iz alignment_info.csv

**5.** Projekcija svih CSV tačaka na pivot poliliniju → stacionaža + signed offset:
```python
def projekcija_na_poliliniju(tacka, polilinija):
    min_dist = np.inf
    best_stacionaza = 0
    best_signed_offset = 0
    kumulativna = 0
    for i in range(len(polilinija) - 1):
        p0 = polilinija[i, :2]
        p1 = polilinija[i+1, :2]
        seg = p1 - p0
        seg_len = np.linalg.norm(seg)
        if seg_len < 1e-10:
            continue
        t = np.dot(tacka[:2] - p0, seg) / seg_len**2
        t = np.clip(t, 0, 1)
        proj = p0 + t * seg
        dist = np.linalg.norm(tacka[:2] - proj)
        if dist < min_dist:
            min_dist = dist
            best_stacionaza = kumulativna + t * seg_len
            tangenta = seg / seg_len
            normala = np.array([-tangenta[1], tangenta[0]])
            best_signed_offset = np.dot(tacka[:2] - proj, normala)
        kumulativna += seg_len
    return best_stacionaza, best_signed_offset
```
Koristiti `np.searchsorted` za efikasno lociranje presjeka (ne KDTree).

**6.** 1D interpolacija Z po presjecima:
- Za svaki presjek: tačke unutar ±2.5m stacionaže
- Sortiraj po signed offsetu → 1D interpolacija Z(offset)
- Geometrijski konzistentno, izbjegava artefakte globalnog 2D interpolatora

**7.** Signed nagibi iz 1D interpolacije: `nagib = (Z_spoljna - Z_pivot) / sirina`

**8.** Smooth ciljnih nagiba po dionicama:
- PRAVAC/KRIVINA: rolling average, prozor = `min(max(duzina_dionice/3, 50m), 150m)`
  - **VAŽNO:** ako je `window > duzina_dionice` → eksplicitno clampuj prozor na dužinu dionice
  - Rastući/padajući trend nagiba koji ostaje iznad `|2.5%|` nije šum — to je geometrija krivine
- PRIJELAZNICA: linearna interpolacija između susjednih dionica (monotona promjena, prolazi kroz nulu)
  - **Edge case:** prijelaznica na rubu koridora (početak/kraj trase) → ekstrapolacija iz jedne susjedne dionice

**9.** Smooth pivot Z — QP sa hard constraints:
- Faza 1: agresivno smooth bez constraint → idealni rezultat
- Faza 2: clip na [-1/+3cm] granice
- Gdje Faza 2 značajno odstupa od Faze 1 → MANUAL CHECK kandidati
- Adaptive prozor: veći gdje je lokalni std_dev mali, manji gdje je std_dev veliki

**10.** Za svaki presjek (na svakih 5m):

**AKO PRAVAC ili KRIVINA:**
- a. Nova pivot Z = smoothed, unutar [-1/+3cm]
- b. Nova spoljna Z = Z_pivot_novo - (nagib_cilj/100 × sirina)
- c. Constraint provjera spoljne Z:
  - `OK`: unutar [-1/+3cm] i `abs(nagib) ≥ 2.5%`
  - `AUTO-ADJUSTED`: van prozora ili `abs(nagib) < 2.5%` → prilagodi nagib prema granici (postavi na `sign × 2.5%`)
  - `MANUAL CHECK`: ni granični nagib ne prolazi
- d. Srednje tačke: linearna interpolacija prema signed offsetu

**AKO PRIJELAZNICA:**
- a. Nova pivot Z = smoothed, unutar [-1/+3cm]
- b. Nagib = linearna interpolacija između nagiba susjednih dionica (monotona, kroz nulu)
- c. Constraint provjera spoljne Z (bez odvodnja checka)
- d. Status = PRIJELAZNICA
- e. Srednje tačke: linearna interpolacija prema signed offsetu

**11.** Output novi_dtm.csv + qa_report.csv

---

## Tehnički stack

```python
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d
import cvxpy as cp   # QP solver za Skriptu 0
# OSQP je backend koji cvxpy koristi automatski
```

---

## Prioriteti implementacije

1. **Skripta 0 — QP pivot smoothing** (temelj za Skriptu 1)
2. **Robustna projekcija sa signed offsetom** (temelj Skripte 1)
3. **Stabilne normale** (reparametrizacija + širi prozor)
4. **Outlier detekcija** (rolling median, bez external alignment info)
5. **Tretman prijelaznica** (monoton nagib, suspenzija odvodnja checka)
6. **1D interpolacija po presjecima** (ne globalni 2D)
7. **Signed nagib kroz cijeli pipeline** (ne apsolutna vrijednost)
8. **Adaptive smoothing + AUTO-ADJUSTED logika**

---

## Cilj kvaliteta

Od ~600 presjeka: svesti MANUAL CHECK na ~20-50.
Finalni vizualni pregled outputa radi se u Civil 3D / ProVi CAD.

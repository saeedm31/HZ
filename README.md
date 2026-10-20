# 🌍 HZ: Habitable Zone Calculator

**Where around a star could a planet have liquid water? Give it a star, get the habitable-zone boundaries in AU.**

Implements the habitable-zone model of [Kopparapu et al. (2013)](https://arxiv.org/abs/1301.6674) for main-sequence stars, from my exoplanet research.

## What it computes

Five boundaries, from the inner to the outer edge:

| Boundary | Meaning |
|---|---|
| Recent Venus | Optimistic inner edge |
| Runaway greenhouse | Inner edge, where oceans would boil away |
| Moist greenhouse | Conservative inner edge |
| Maximum greenhouse | Conservative outer edge |
| Early Mars | Optimistic outer edge |

These are combined into a **conservative** range (moist to maximum greenhouse) and an **optimistic** one (recent Venus to early Mars).

## Usage

**From effective temperature and luminosity** (`HZ_TeffLum_Main.py`):

```python
from HZ_TeffLum_Main import HZ_in_out_2013

rv, rg, mg, maxg, em, conservative, optimistic = HZ_in_out_2013(T_eff=5780, Lum=1.0)
print(conservative, optimistic)   # Sun: [0.99, 1.71] AU and [0.75, 1.77] AU
```

**From effective temperature only** (`HZ_MassTeff.py`): when the luminosity is unknown, it is estimated from a stellar model grid (`stars-model.txt`). Works for T_eff ≥ 3050 K.

```python
from HZ_MassTeff import HZ_InOut_MassTeff

*_, conservative, optimistic, lum = HZ_InOut_MassTeff(0.8, 5000)
print(conservative, lum)          # 5000 K star: [0.57, 1.01] AU, L ≈ 0.31 L☉
```

Temperatures in kelvin, luminosity in solar units, distances in AU. Requires `numpy` and `scipy`.

## Example results

| Star | T_eff | L (L☉) | Conservative HZ (AU) | Optimistic HZ (AU) |
|---|---|---|---|---|
| Sun | 5780 K | 1.0 | 0.99 – 1.71 | 0.75 – 1.77 |
| K dwarf | 5000 K | 0.31 (estimated) | 0.57 – 1.01 | 0.43 – 1.05 |
| M dwarf | 3300 K | 0.0125 | 0.12 – 0.23 | 0.09 – 0.24 |

---

Part of my astronomy background: [Saeed Hojjatpanah](https://github.com/saeedm31), PhD in Astronomy (exoplanets).

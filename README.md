# tafilah-wind-farm-design
90 MW wind farm design (Al-Tafilah, Jordan) - crescent-arc layout for 30 Vestas V112 turbines. Compares Jensen (Park) and Gaussian wake models in Python. AEP: 256- 283 GWh/yr, wake loss: 11.6-20.1%.
# 90 MW Wind Farm Design – Al-Tafilah, Jordan

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Design and simulation of a 90 MW wind farm using 30 Vestas V112 turbines in a crescent-arc layout. The project compares two analytical wake models:

- **Jensen (Park)** model (top‑hat profile, worst‑case superposition)
- **Gaussian wake model** (Bastankhah & Porté-Agel, 2014) with root‑sum‑square superposition

Annual energy production (AEP) is computed using 8,760 hours of synthetic wind data derived from monthly averages at the site.

## Key Results

| Wake Model       | AEP (GWh/yr) | Wake Loss (%) | Capacity Factor |
|-----------------|--------------|---------------|------------------|
| No Wake (baseline) | 320.30       | –             | 40.6%            |
| Jensen (Park)   | 283.28       | 11.56%        | 35.9%            |
| Gaussian        | 256.09       | 20.05%        | 32.5%            |

The crescent layout effectively sheds wakes from outer‑ring turbines; inner‑ring turbines experience up to 33.7% wake loss (Gaussian model).

## Repository Contents

- `wind_farm_sim.py` – Python simulation code (Jensen and Gaussian models)
- `WindFarmReport.pdf` – Full project report (62 pages)
- `requirements.txt` – Python dependencies
- `LICENSE` – MIT License

## How to Run

```bash
pip install -r requirements.txt
python wind_farm_sim.py

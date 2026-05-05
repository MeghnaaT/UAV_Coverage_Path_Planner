# 🚁 UAV Autonomous Coverage Path Planner

<p align="left">
  <img src="https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white"/>
  <img src="https://img.shields.io/badge/Shapely-Geospatial-green?style=flat"/>
  <img src="https://img.shields.io/badge/Algorithm-A*%20%7C%20TSP%20%7C%20Boustrophedon-blue?style=flat"/>
  <img src="https://img.shields.io/badge/Output-KML%20%7C%20Matplotlib-orange?style=flat"/>
  <img src="https://img.shields.io/badge/License-MIT-lightgrey?style=flat"/>
</p>

An autonomous survey coverage system for UAVs that computes optimal lawnmower flight paths over arbitrary polygon regions. Given a KML boundary file, the system decomposes the area into convex sub-regions, generates parallel sweep lines calibrated to camera FOV and overlap, and stitches them into a complete mission path using A\* transitions and TSP ordering — minimizing total flight distance while guaranteeing area coverage.

---

## 📽️ Demo Output

> *Matplotlib visualization and KML output — importable directly into Google Earth / mission planners.*

The system outputs two KML files:
- `optimized_convex_decomposition.kml` — color-coded convex sub-regions overlaid on input boundary
- `lawnmower_paths.kml` — complete stitched mission path with inter-cell transitions

---

## ✨ Key Features

- **Arbitrary polygon support** — handles non-convex, irregular survey boundaries via KML input
- **Hertel-Mehlhorn convex decomposition** — decomposes non-convex polygons into minimal convex sub-cells ensuring complete, overlap-free coverage
- **Camera-aware path spacing** — automatically computes sweep line spacing from drone altitude, lens FOV (horizontal + vertical), camera pitch, and desired overlap percentage
- **Boustrophedon (lawnmower) sweep** — generates parallel sweep lines anchored to the longest edge of each convex cell for efficient coverage
- **A\* inter-cell transitions** — computes obstacle-aware transition paths between cells constrained to stay within the polygon boundary
- **TSP cell ordering** — solves the Travelling Salesman Problem across cells to minimize total mission distance
- **KML export** — outputs mission paths directly importable into Google Earth Pro and ArduPilot/QGroundControl-compatible planners

---

## 🏗️ System Architecture

```
Input KML (polygon boundary)
        │
        ▼
┌─────────────────────────┐
│  Polygon Extraction     │  minidom KML parser → Shapely Polygon
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Convex Decomposition   │  Hertel-Mehlhorn → merge triangles → convex cells
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Path Spacing Calc      │  Altitude × tan(FOV/2) × (1 - overlap)
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  Lawnmower Path Gen     │  Sweep lines anchored to longest edge, zigzag order
└────────────┬────────────┘
             │
             ▼
┌─────────────────────────┐
│  TSP + A* Stitching     │  Optimal cell ordering + constrained transitions
└────────────┬────────────┘
             │
             ▼
     KML Output + Matplotlib Visualization
```

---

## ⚙️ Camera Configuration

The path spacing is automatically computed from your drone's camera specs:

```python
spacing = get_path_spacing(
    altitude=12,          # meters AGL
    horizontal_fov=130,   # degrees
    vertical_fov=130,     # degrees
    overlap=0.3           # 30% overlap between passes
)
```

> Tested with **Hawkeye Thumb 4K** specs (168.5° H-FOV, 159.76° V-FOV). Adapt constants in `get_path_spacing()` for your sensor.

---

## 🚀 Getting Started

### Prerequisites

```bash
pip install shapely simplekml matplotlib numpy
```

### Usage

1. **Prepare your input KML** — export a polygon boundary from Google Earth Pro as a `.kml` file.

2. **Configure paths and camera in `Path_planning_code.py`:**

```python
# In the __main__ block:
input_kml_path  = "path/to/your/input_polygon.kml"
output_kml_path = "optimized_convex_decomposition.kml"
path_kml_output = "lawnmower_paths.kml"
start_pos       = (longitude, latitude)   # UAV home position

spacing = get_path_spacing(
    altitude=12,
    horizontal_fov=130,
    vertical_fov=130,
    overlap=0.3
)
```

3. **Run:**

```bash
python Path_planning_code.py
```

4. **View outputs** — open generated `.kml` files in Google Earth Pro or import into your mission planner.

---

## 📁 Repository Structure

```
├── Path_planning_code.py       # Main pipeline (all modules)
├── input_polygon.kml           # Sample input boundary (Bhopal test area)
├── lawnmower_paths.kml         # Sample output: stitched mission path
├── optimized_convex_decomposition.kml  # Sample output: decomposed cells
└── README.md
```

---

## 🧮 Algorithms Used

| Component | Algorithm | Why |
|---|---|---|
| Area decomposition | Hertel-Mehlhorn (triangle merging) | Minimal convex partition of arbitrary polygons |
| Sweep direction | Boustrophedon (longest-edge aligned) | Minimizes turn count and path length |
| Cell ordering | TSP via permutation (brute-force for small n) | Optimal inter-cell travel order |
| Transition paths | A\* search constrained to polygon | Avoids leaving survey boundary |
| Coverage metric | Shapely area union ratio | Verifies ≥99% coverage before output |

---

## 📈 Performance Notes

On a sample test polygon (Bhopal region, ~0.001 sq degrees):
- Convex decomposition: typically 4–8 cells for irregular polygons
- Coverage ratio: ≥95% area coverage (verified via Shapely union)
- TSP ordering reduces inter-cell transit vs. naive sequential ordering

---

## 🔮 Planned Improvements

- [ ] Add no-fly zone obstacle support
- [ ] ROS integration for direct autopilot command

---

## 🤝 Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss the proposed change.

1. Fork the repository
2. Create a feature branch: `git checkout -b feat/your-feature`
3. Commit with a descriptive message: `git commit -m "feat: add no-fly zone obstacle avoidance"`
4. Open a PR against `main`

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

*Built by [Meghna Tiwari](https://github.com/MeghnaaT) · Open to feedback and collaboration*


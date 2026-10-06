# Lithography Lessons

Short interactive lessons on the Mo/Si multilayer mirrors used in EUV lithography.

Live site: https://photolitho.github.io/lithography-lessons/

| Module | Question | Status |
|---|---|---|
| 1 | How much EUV light does a single surface reflect? | Available |
| 2 | What do two surfaces close together do? | Available |
| 3 | How do 80 weak surfaces make a 74% mirror? | Available |
| 4 | How narrow is the peak, and what moves it? | Available |
| 5 | Where does the light that does not come back go? | Available |
| 6 | How hot does the mirror get? | Available |
| 7 | How much does the hot surface move? | Planned |
| 8 | How should the mirror be cooled? | Planned |

The home page opens on a cover; after a touch the mirror builds itself in five seconds and becomes the map of the modules: each band of the drawing opens its module and shows the number it produces. Below the map, the page explains the project and shows the chain of mirrors and the mirror in two drawings. The opening sequence uses GSAP 3.12.5 from cdnjs; without it the page shows the map directly. `manifest.webmanifest` declares the project folder as one home-screen web app.

## How it is built

Each module is one HTML file with its own text, Italian dictionary and experiments. Shared by all modules: `lesson.css`, `lesson.js` (navigation, language switch, presets, drawing and physics helpers) and `cxro.js` (the optical-constant tables). Module 6 also loads `thermal.js` (Bessel functions and the steady temperature of a heated disc). Charts are drawn on canvas. There is no build step and no dependency beyond web fonts.

Every module has a single view with nothing to answer or complete (`open: true` in its `startLesson` call). Its first page is an index of its ideas. Every idea opens with a short "Why here" paragraph that says what it builds on and what uses it later, then shows its experiment, its explanation and a card for each important formula: where it comes from, its assumptions, the meaning of its terms and one check the reader can do. Experiments offer a few presets ("Try") that set a telling case and say what to look at. Mathematical prerequisites are introduced at the idea that first needs them. The last visited idea is kept in the browser's local storage.

## Languages

Lessons are available in English and Italian. The EN / IT switch at the top changes language on the spot, keeps the current step and answers, and is remembered across pages. English is the source: Italian strings live in a dictionary keyed by the English text, and anything missing falls back to English.

## Data

Optical constants of Mo and Si come from the CXRO database (Henke, Gullikson and Davis, Lawrence Berkeley National Laboratory). The tables are embedded exactly as tabulated (photon energy, delta, beta) and interpolated linearly in wavelength, with `n = 1 - delta` and `k = beta`.

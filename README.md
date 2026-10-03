# Lithography Lessons

Short interactive lessons on the Mo/Si multilayer mirrors used in EUV lithography.

Live site: https://photolitho.github.io/lithography-lessons/

| Module | Question | Status |
|---|---|---|
| 1 | How much EUV light does a single surface reflect? | Available |
| 2 | What do two surfaces close together do? | Available |
| 3 | How do 80 weak surfaces make a 74% mirror? | Available |
| 4 | Wavelength, period and the share of molybdenum | In preparation |
| 5 | Where the light that does not come back goes | In preparation |

## How it is built

Each module is one HTML file with its own text, Italian dictionary and experiments. Three files are shared by all modules: `lesson.css`, `lesson.js` (navigation, language switch, questions, drawing and physics helpers) and `cxro.js` (the optical-constant tables). Charts are drawn on canvas. There is no build step and no dependency beyond web fonts.

Every important formula comes with a card that states where it comes from, its assumptions, the meaning of its terms and one check the reader can do. Mathematical prerequisites are introduced at the step that first needs them.

Every module has two modes. The guided path asks for a prediction and unlocks each step with a question or a goal. Explore gives direct access to all steps with no required answers. Progress is kept separately for the two modes, in the browser's local storage.

## Languages

Lessons are available in English and Italian. The EN / IT switch at the top changes language on the spot, keeps the current step and answers, and is remembered across pages. English is the source: Italian strings live in a dictionary keyed by the English text, and anything missing falls back to English.

## Data

Optical constants of Mo and Si come from the CXRO database (Henke, Gullikson and Davis, Lawrence Berkeley National Laboratory). The tables are embedded exactly as tabulated (photon energy, delta, beta) and interpolated linearly in wavelength, with `n = 1 - delta` and `k = beta`.

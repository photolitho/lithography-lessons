# Lithography Lessons

Short interactive lessons on the Mo/Si multilayer mirrors used in EUV lithography.

Live site: https://photolitho.github.io/lithography-lessons/

| Module | Question | Status |
|---|---|---|
| 1 | How much EUV light does a single surface reflect? | Available |
| 2 | A single layer | In preparation |
| 3 | The full stack | In preparation |

## How it is built

Each module is one self-contained HTML file: markup, styles, script and data together, with charts drawn on canvas. There is no build step and no dependency beyond web fonts.

Every module has two modes. The guided path asks for a prediction and unlocks each step with a question or a goal. Explore gives direct access to all steps with no required answers. Progress is kept separately for the two modes, in the browser's local storage.

## Languages

Lessons are available in English and Italian. The EN / IT switch at the top changes language on the spot, keeps the current step and answers, and is remembered across pages. English is the source: Italian strings live in a dictionary keyed by the English text, and anything missing falls back to English.

## Data

Optical constants of Mo and Si come from the CXRO database (Henke, Gullikson and Davis, Lawrence Berkeley National Laboratory). The tables are embedded exactly as tabulated (photon energy, delta, beta) and interpolated linearly in wavelength, with `n = 1 - delta` and `k = beta`.

// Shared engine for the lithography lessons: physics helpers, language, drawing helpers, navigation.
const LAM0 = 13.5;

/* ---------- physics ---------- */
const C = {
  of: x => Array.isArray(x) ? x : [x, 0],
  add: (a, b) => [a[0] + b[0], a[1] + b[1]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
  sqrt: a => { const m = Math.hypot(a[0], a[1]); const re = Math.sqrt((m + a[0]) / 2); const im = Math.sqrt(Math.max(0, (m - a[0]) / 2)); return [re, a[1] < 0 ? -im : im]; },
  abs2: a => a[0] * a[0] + a[1] * a[1]
};
function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let lo = 0, hi = xs.length - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] > x) hi = m; else lo = m; }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}
const index = (mat, lam) => [interp(DATA[mat].l, DATA[mat].n, lam), interp(DATA[mat].l, DATA[mat].k, lam)];
// Amplitude reflection coefficient; theta = angle from the normal in medium 1.
function fresnel(n1, n2, theta, pol) {
  n1 = C.of(n1); n2 = C.of(n2);
  const c1 = [Math.cos(theta), 0], s2 = C.mul(C.div(n1, n2), [Math.sin(theta), 0]);
  const c2 = C.sqrt(C.sub([1, 0], C.mul(s2, s2)));
  const [a, b] = pol === "s" ? [C.mul(n1, c1), C.mul(n2, c2)] : [C.mul(n2, c1), C.mul(n1, c2)];
  return C.div(C.sub(a, b), C.add(a, b));
}
const attLen = (k, lam) => lam / (4 * Math.PI * k);
// ----- layered stacks (transfer-matrix method) -----
// q = n·cos θ inside a material, from Snell's law with the vacuum-side angle; root chosen so that waves decay with depth.
function qOf(n, n0sin) { n = C.of(n); const q = C.sqrt(C.sub(C.mul(n, n), [n0sin * n0sin, 0])); return q[1] < 0 ? [-q[0], -q[1]] : q; }
// Admittance: the ratio H/E of a forward wave. s: q, p: n²/q.
const admittance = (n, q, pol) => pol === "s" ? q : C.div(C.mul(C.of(n), C.of(n)), q);
const cCos = z => [Math.cos(z[0]) * Math.cosh(z[1]), -Math.sin(z[0]) * Math.sinh(z[1])];
const cSin = z => [Math.sin(z[0]) * Math.cosh(z[1]), Math.cos(z[0]) * Math.sinh(z[1])];
const mMul = (A, B) => [C.add(C.mul(A[0], B[0]), C.mul(A[1], B[2])), C.add(C.mul(A[0], B[1]), C.mul(A[1], B[3])), C.add(C.mul(A[2], B[0]), C.mul(A[3], B[2])), C.add(C.mul(A[2], B[1]), C.mul(A[3], B[3]))];
// Matrix of one layer, linking (E, H) at its top to (E, H) at its bottom. Stored as [m11, m12, m21, m22].
function layerMatrix(n, d, lam, n0sin, pol) {
  const q = qOf(n, n0sin), eta = admittance(n, q, pol), delta = C.mul(q, [2 * Math.PI * d / lam, 0]), c = cCos(delta), s = cSin(delta), mi = [0, -1];
  return [c, C.div(C.mul(mi, s), eta), C.mul(C.mul(mi, eta), s), c];
}
// Reflection and transmission of a stack. layers: [{ n, d }] from the vacuum side down; n0, ns: indices above and below.
// Returns r (amplitude), R, T (power entering the substrate) and A = 1 − R − T (power absorbed in the layers).
function tmm(layers, n0, ns, lam, theta0, pol) {
  const n0sin = C.of(n0)[0] * Math.sin(theta0), eta0 = admittance(n0, qOf(n0, n0sin), pol), etaS = admittance(ns, qOf(ns, n0sin), pol);
  let M = [[1, 0], [0, 0], [0, 0], [1, 0]];
  for (const L of layers) M = mMul(M, layerMatrix(L.n, L.d, lam, n0sin, pol));
  const B = C.add(M[0], C.mul(M[1], etaS)), Cc = C.add(M[2], C.mul(M[3], etaS)), den = C.add(C.mul(eta0, B), Cc);
  const r = C.div(C.sub(C.mul(eta0, B), Cc), den), tt = C.div(C.mul([2, 0], eta0), den), R = C.abs2(r), T = etaS[0] / eta0[0] * C.abs2(tt);
  return { r, R, T, A: 1 - R - T, M, eta0, etaS, n0sin };
}
// Fields and energy flow inside a stack, for an incoming wave of amplitude 1.
// Returns, besides r, R, T: tops[j] = (E, H) at the top of layer j (last entry: top of the substrate),
// flux[j] = net power crossing that plane as a fraction of the incident power, absorbed[j] = flux[j] − flux[j+1].
function stackFields(layers, n0, ns, lam, theta0, pol) {
  const x = tmm(layers, n0, ns, lam, theta0, pol), eta0 = x.eta0;
  let E = C.add([1, 0], x.r), H = C.mul(eta0, C.sub([1, 0], x.r));
  const info = layers.map(L => { const q = qOf(L.n, x.n0sin); return { q, eta: admittance(L.n, q, pol), d: L.d }; });
  const fluxOf = (e, h) => (e[0] * h[0] + e[1] * h[1]) / eta0[0];      // Re(E·conj(H)) / Re(η0)
  // State a depth z below a plane, inside a material: the inverse of the layer matrix.
  const down = (e, h, L, z) => { const dl = C.mul(L.q, [2 * Math.PI * z / lam, 0]), c = cCos(dl), sn = cSin(dl), i = [0, 1];
    return [C.add(C.mul(c, e), C.div(C.mul(C.mul(i, sn), h), L.eta)), C.add(C.mul(C.mul(C.mul(i, L.eta), sn), e), C.mul(c, h))]; };
  const tops = [], flux = [];
  for (const L of info) { tops.push([E, H]); flux.push(fluxOf(E, H)); [E, H] = down(E, H, L, L.d); }
  tops.push([E, H]); flux.push(fluxOf(E, H));
  const absorbed = info.map((L, j) => flux[j] - flux[j + 1]);
  // |E|² at depth z below the top of layer j.
  const E2 = (j, z) => C.abs2(down(tops[j][0], tops[j][1], info[j], z)[0]);
  return Object.assign(x, { tops, flux, absorbed, E2, info });
}
// N pairs of two materials, top layer first.
const pairs = (N, a, b) => Array.from({ length: 2 * N }, (_, i) => i % 2 ? b : a);
// Factor a wave picks up crossing a film of complex index n and thickness d twice: e^(i·4π·n·d/λ).
function roundTrip(n, d, lam) { const a = 4 * Math.PI * d / lam, m = Math.exp(-a * n[1]); return [m * Math.cos(a * n[0]), m * Math.sin(a * n[0])]; }

/* ---------- language ---------- */
// Italian strings keyed by their English source; anything missing falls back to English.
// Strings shared by every module live here; each module adds its own with Object.assign(IT, {...}).
const IT = {
  "Menu": "Menu",
  "← All lessons": "← Tutte le lezioni",
  "Guided path": "Percorso guidato",
  "One opening prediction, then one step at a time. You move on by answering a question or reaching a goal.": "Una previsione iniziale, poi un passo alla volta. Si avanza rispondendo a una domanda o raggiungendo un obiettivo.",
  "Explore": "Esplora",
  "Every experiment, in any order, with no required answers. Explanations open on request.": "Tutti gli esperimenti, in qualsiasi ordine, senza risposte obbligatorie. Le spiegazioni si aprono su richiesta.",
  "Exploring does not complete steps of the guided path.": "Esplorare non completa i passi del percorso guidato.",
  "Back": "Indietro",
  "Check": "Verifica",
  "Step": "Passo",
  "Hide explanation": "Nascondi la spiegazione",
  "Show explanation": "Mostra la spiegazione",
  "Continue": "Continua",
  "Module complete": "Modulo completato",
  "Last step": "Ultimo passo",
  "{n} things left to do": "Restano {n} cose da fare",
  "One thing left to do": "Resta una cosa da fare",
  "no prediction": "nessuna previsione",
  "{n} of {m} steps completed": "{n} di {m} passi completati",
  "Start": "Inizia",
  "Review": "Rivedi",
  "Resume": "Riprendi",
  "Correct.": "Esatto.",
  "Type a number.": "Scrivi un numero.",
  "Where it comes from": "Da dove viene",
  "Assumptions": "Ipotesi",
  "The terms": "I termini",
  "vacuum": "vuoto",
  "material": "materiale",
  "silicon": "silicio",
  "molybdenum": "molibdeno",
  "Silicon": "Silicio",
  "Molybdenum": "Molibdeno",
  "Module": "Modulo",
  "Guided": "Guidato",
  "EUV Mirror": "Specchio EUV",
  "EUV Module": "EUV Modulo"
};
const LANG_KEY = "litho-lang";
let lang = "en";
try { lang = localStorage.getItem(LANG_KEY) || ((navigator.language || "").toLowerCase().startsWith("it") ? "it" : "en"); } catch (e) {}
const t = s => (lang === "it" && IT[s]) || s;
const loc = v => lang === "it" ? String(v).replace(".", ",") : String(v);

/* ---------- helpers ---------- */
const $ = id => document.getElementById(id);
const fmt = (x, d) => loc(x.toFixed(d));
const pct = (x, d) => fmt(100 * x, d) + "%";
// Percentage with digits suited to the magnitude, for values spanning many orders.
const pctAuto = x => { const p = 100 * x; return (p >= 10 ? fmt(p, 1) : p >= 0.01 ? fmt(p, 3) : loc(Number(p.toPrecision(2)))) + "%"; };
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function prep(c) {
  const r = c.getBoundingClientRect(), d = window.devicePixelRatio || 1;
  c.width = Math.round(r.width * d); c.height = Math.round(r.height * d);
  const ctx = c.getContext("2d"); ctx.setTransform(d, 0, 0, d, 0, 0); ctx.clearRect(0, 0, r.width, r.height);
  return { ctx, w: r.width, h: r.height };
}
const PAD = { l: 60, r: 12, t: 10, b: 44 };
function plot(c, o) {
  const { ctx, w, h } = prep(c);
  const ty = v => o.logy ? Math.log10(v) : v;
  const X = x => PAD.l + (x - o.x[0]) / (o.x[1] - o.x[0]) * (w - PAD.l - PAD.r);
  const Y = y => h - PAD.b - (ty(y) - ty(o.y[0])) / (ty(o.y[1]) - ty(o.y[0])) * (h - PAD.t - PAD.b);
  ctx.font = "11px " + css("--mono"); ctx.fillStyle = css("--muted"); ctx.strokeStyle = css("--line"); ctx.lineWidth = 1;
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of o.xt) { ctx.beginPath(); ctx.moveTo(X(t), PAD.t); ctx.lineTo(X(t), h - PAD.b); ctx.stroke(); ctx.fillText(o.xf ? o.xf(t) : t, X(t), h - PAD.b + 6); }
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const t of o.yt) { ctx.beginPath(); ctx.moveTo(PAD.l, Y(t)); ctx.lineTo(w - PAD.r, Y(t)); ctx.stroke(); ctx.fillText(o.yf ? o.yf(t) : t, PAD.l - 7, Y(t)); }
  ctx.save(); ctx.beginPath(); ctx.rect(PAD.l, PAD.t, w - PAD.l - PAD.r, h - PAD.t - PAD.b); ctx.clip();
  for (const bd of o.bands || []) { ctx.globalAlpha = .22; ctx.fillStyle = bd.color; ctx.fillRect(X(bd.x0), PAD.t, X(bd.x1) - X(bd.x0), h - PAD.t - PAD.b); ctx.globalAlpha = 1; }
  if (o.vline != null) { ctx.strokeStyle = css("--accent"); ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(X(o.vline), PAD.t); ctx.lineTo(X(o.vline), h - PAD.b); ctx.stroke(); ctx.setLineDash([]); }
  for (const s of o.series) {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.width || 2.2; ctx.setLineDash(s.dash || []);
    ctx.beginPath(); s.pts.forEach((p, i) => i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1]))); ctx.stroke(); ctx.setLineDash([]);
  }
  for (const m of o.marks || []) { ctx.fillStyle = m.color; ctx.beginPath(); ctx.arc(X(m.x), Y(m.y), 5, 0, 7); ctx.fill(); }
  ctx.restore();
  if (o.xlabel) { ctx.fillStyle = css("--muted"); ctx.textAlign = "center"; ctx.textBaseline = "bottom"; ctx.fillText(o.xlabel, (PAD.l + w - PAD.r) / 2, h - 2); }
}
// Arrow from (x0, y0) to (x1, y1), used for phasors and rays.
function arrow(ctx, x0, y0, x1, y1, color, width) {
  const len = Math.hypot(x1 - x0, y1 - y0), a = Math.atan2(y1 - y0, x1 - x0), hd = Math.min(11, len * .6);
  ctx.strokeStyle = ctx.fillStyle = color; ctx.lineWidth = width; ctx.lineCap = "round";
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
  if (len < 3) return;
  ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x1 - hd * Math.cos(a - .42), y1 - hd * Math.sin(a - .42)); ctx.lineTo(x1 - hd * Math.cos(a + .42), y1 - hd * Math.sin(a + .42)); ctx.closePath(); ctx.fill();
}
// Strip of material: light enters from the left and fades.
function strip(c, frac, color, inset) {
  const { ctx, w, h } = prep(c), x0 = inset ? PAD.l : 0, x1 = w - (inset ? PAD.r : 0);
  for (let x = x0; x < x1; x++) { ctx.globalAlpha = Math.max(0, frac((x - x0) / (x1 - x0))); ctx.fillStyle = color; ctx.fillRect(x, 0, 1.2, h); }
  ctx.globalAlpha = 1; ctx.strokeStyle = css("--line"); ctx.strokeRect(x0 + .5, .5, x1 - x0 - 1, h - 1);
}

/* ---------- lesson engine ---------- */
// Shared by every module. A module page supplies its markup, its Italian strings (Object.assign(IT, {...}))
// and its experiments (the setup function), then calls startLesson({ n, key, next, setup }).
// Two modes with separate progress: "guided" uses cur and solved, "explore" uses ecur and never writes to solved.
let cfg, app, home, steps, st, tEls;
const draws = {}, opened = {};
// Static text to translate: outermost matches only, with their English HTML kept as the source.
const TSEL = "h1, h2, p, summary, dt, dd, .carry li, .cap:not([id]), .goal, .ctrl label, .opts button, .toggle button, .reads b, [data-t], #back, .numq .btn, #menu, a.linkbtn";
const save = () => { try { localStorage.setItem(cfg.key, JSON.stringify(st)); } catch (e) {} };
const setFb = (fb, en) => { fb.dataset.en = en; fb.textContent = t(en); };
const guided = () => st.mode === "guided";
// Steps are tracked by their data-id, so adding or reordering steps does not scramble saved progress.
const idOf = i => steps[i].dataset.id;
const stepOf = el => steps.indexOf(el.closest(".step"));
const isSolved = i => st.solved.includes(idOf(i));
const reach = () => { let i = 0; while (i < steps.length - 1 && isSolved(i)) i++; return i; };
const pos = () => Math.min(guided() ? st.cur : st.ecur, steps.length - 1);
const title = i => steps[i].querySelector("h1").textContent;
const stepLabel = i => t("Step") + " " + (i + 1) + ": " + title(i);

// The explanation shows once the step is solved (guided) or on request (explore).
function syncAfter(i) {
  const a = steps[i].querySelector(".after"); if (!a) return;
  a.hidden = guided() ? !isSolved(i) : !opened[i];
  steps[i].querySelector(".expl").textContent = t(opened[i] ? "Hide explanation" : "Show explanation");
}
// A step is solved when every question and goal in it is done.
const reqs = i => [...steps[i].querySelectorAll(".q, .goal")];
function done(el) {
  if (!guided()) return;
  el.classList.add("done");
  const i = stepOf(el);
  if (!isSolved(i) && reqs(i).every(r => r.classList.contains("done"))) { st.solved.push(idOf(i)); save(); }
  syncAfter(i); refresh();
}
function refresh() {
  const g = guided(), cur = pos(), r = reach(), last = cur === steps.length - 1;
  [...$("segs").children].forEach((b, i) => { b.className = "seg" + (g && isSolved(i) ? " done" : "") + (i === cur ? " cur" : "") + (!g || i <= r ? " reach" : ""); });
  $("where").textContent = t("Module") + " " + cfg.n + " · " + t(g ? "Guided" : "Explore");
  $("count").textContent = (cur + 1) + " / " + steps.length;
  $("back").disabled = cur === 0;
  $("next").disabled = last || (g && !isSolved(cur));
  $("next").textContent = t(!last ? "Continue" : g ? "Module complete" : "Last step");
  const left = reqs(cur).filter(r => !r.classList.contains("done")).length;
  $("hint").textContent = !g ? "" : !isSolved(cur) ? t(left > 1 ? "{n} things left to do" : "One thing left to do").replace("{n}", left) : last ? t(cfg.next || "") : "";
}
function show(i, keep) {
  if (guided()) st.cur = i; else st.ecur = i;
  save();
  home.hidden = true; steps.forEach((s, j) => s.hidden = j !== i);
  if ($("pred-echo")) $("pred-echo").textContent = t(st.pred || "no prediction");
  syncAfter(i); refresh(); if (draws[idOf(i)]) draws[idOf(i)]();
  if (!keep) window.scrollTo(0, 0);
}
// Put sliders and toggles back to their defaults, so a goal is never met by what was done in the other mode.
function resetControls() {
  document.querySelectorAll("input[type=range]").forEach(r => r.value = r.defaultValue);
  document.querySelectorAll(".toggle").forEach(tg => tg.firstElementChild.click());
}
function enter(mode, i) {
  st.mode = mode; app.dataset.mode = mode; resetControls();
  show(i != null ? i : mode === "guided" ? Math.min(st.cur, reach()) : pos());
}
function goHome(keep) {
  st.mode = "home"; app.dataset.mode = "home"; save();
  steps.forEach(s => s.hidden = true); home.hidden = false;
  const n = steps.filter((s, i) => isSolved(i)).length;
  $("home-prog").textContent = t("{n} of {m} steps completed").replace("{n}", n).replace("{m}", steps.length);
  $("go-guided").textContent = t(n === 0 && !st.pred ? "Start" : n === steps.length ? "Review" : "Resume");
  $("where").textContent = t("EUV Mirror") + " · " + t("Module") + " " + cfg.n;
  if (!keep) window.scrollTo(0, 0);
}
// Experiments call goal() on every redraw; it only counts in guided mode, once.
function goal(id, ok) { const g = $(id); if (ok && guided() && !g.classList.contains("done")) done(g); }
function toggle(id, cb) { const bs = [...$(id).children]; bs.forEach(b => b.onclick = () => { bs.forEach(x => x.setAttribute("aria-pressed", x === b)); cb(b.dataset); }); }
const redraw = () => { const d = st.mode !== "home" && draws[idOf(pos())]; if (d) d(); };
// Rewrites every static text in the current language; dynamic text follows on the next render.
function translateDom() {
  document.documentElement.lang = lang; document.title = t("EUV Module") + " " + cfg.n;
  tEls.forEach(([e, en]) => e.innerHTML = t(en));
  document.querySelectorAll(".fb").forEach(f => { if (f.dataset.en) f.textContent = t(f.dataset.en); });
  steps.forEach((s, i) => { $("toc").children[i].firstChild.textContent = title(i); $("segs").children[i].setAttribute("aria-label", stepLabel(i)); });
  [...$("t-lang").children].forEach(b => b.setAttribute("aria-pressed", b.dataset.l === lang));
}

function startLesson(config) {
  cfg = config; app = document.querySelector(".app"); home = $("home");
  steps = [...document.querySelectorAll(".step")].filter(s => s !== home);
  st = { mode: "home", cur: 0, ecur: 0, solved: [], pred: null };
  try { Object.assign(st, JSON.parse(localStorage.getItem(cfg.key) || "{}")); } catch (e) {}
  tEls = [...document.querySelectorAll(TSEL)].filter(e => !e.parentElement.closest(TSEL)).map(e => [e, e.innerHTML]);

  steps.forEach((s, i) => {
    const b = document.createElement("button"); b.className = "seg"; b.setAttribute("aria-label", stepLabel(i));
    b.onclick = () => { if (!guided() || i <= reach()) show(i); }; $("segs").append(b);
    const li = document.createElement("li"), tb = document.createElement("button"); tb.textContent = title(i); tb.onclick = () => enter("explore", i); li.append(tb); $("toc").append(li);
    const a = s.querySelector(".after");
    if (a) { const e = document.createElement("button"); e.className = "btn expl explore-only"; e.onclick = () => { opened[i] = !opened[i]; syncAfter(i); }; a.before(e); }
  });
  $("go-guided").onclick = () => enter("guided");
  $("menu").onclick = () => goHome();
  $("back").onclick = () => show(pos() - 1);
  $("next").onclick = () => show(pos() + 1);

  // Multiple-choice questions
  document.querySelectorAll(".mc").forEach(mc => {
    const q = mc.closest(".q"), fb = q.querySelector(".fb"), btns = [...mc.children];
    const win = () => { btns.forEach(b => { b.disabled = true; if ("ok" in b.dataset) b.classList.add("right"); }); setFb(fb, "Correct."); fb.classList.add("ok"); };
    btns.forEach(b => b.onclick = () => { if ("ok" in b.dataset) { win(); done(q); } else { b.classList.add("wrong"); b.disabled = true; setFb(fb, b.dataset.hint); } });
    if (isSolved(stepOf(q))) { win(); q.classList.add("done"); }
  });
  // Numeric questions: data-ans is the expected value, data-tol the accepted distance from it.
  document.querySelectorAll(".numq").forEach(nq => {
    const q = nq.closest(".q"), fb = q.querySelector(".fb"), inp = nq.querySelector("input"), btn = nq.querySelector("button");
    const ans = +nq.dataset.ans, tol = +nq.dataset.tol;
    const win = () => { setFb(fb, "Correct."); fb.classList.add("ok"); };
    btn.onclick = () => {
      const v = parseFloat(inp.value.replace(",", ".").replace("−", "-").replace("%", ""));
      fb.classList.remove("ok");
      if (isNaN(v)) setFb(fb, "Type a number.");
      else if (Math.abs(v - ans) <= tol) { win(); done(q); }
      else setFb(fb, nq.dataset.hint);
    };
    inp.onkeydown = e => { if (e.key === "Enter") btn.click(); };
    if (isSolved(stepOf(q))) { inp.value = loc(ans); win(); q.classList.add("done"); }
  });
  steps.forEach((s, i) => { if (isSolved(i)) s.querySelectorAll(".goal").forEach(g => g.classList.add("done")); });
  // Opening prediction: any choice counts, and it is echoed later in the element with id pred-echo.
  if ($("pred")) [...$("pred").children].forEach(b => {
    const q = $("pred").closest(".q");
    if (st.pred === b.dataset.v) { b.classList.add("picked"); if (isSolved(stepOf(q))) q.classList.add("done"); }
    b.onclick = () => { [...$("pred").children].forEach(x => x.classList.remove("picked")); b.classList.add("picked"); st.pred = b.dataset.v; done(q); };
  });

  if (cfg.setup) cfg.setup();

  window.addEventListener("resize", redraw);
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(redraw);
  [...$("t-lang").children].forEach(b => b.onclick = () => {
    lang = b.dataset.l; try { localStorage.setItem(LANG_KEY, lang); } catch (e) {}
    translateDom(); if (st.mode === "home") goHome(true); else show(pos(), true);
  });
  translateDom();
  if (st.mode === "guided" || st.mode === "explore") enter(st.mode); else goHome();
}

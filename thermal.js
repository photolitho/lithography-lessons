// Heat conduction helpers for the thermal modules.
// Bessel functions J0, J1 (rational approximations, Numerical Recipes), accurate to about 1e-8.
function besselJ0(x) { const ax = Math.abs(x);
  if (ax < 8) { const y = x * x; return (57568490574 + y * (-13362590354 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * -184.9052456))))) / (57568490411 + y * (1029532985 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 0.785398164;
  return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))) - z * Math.sin(xx) * (-0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7))))); }
function besselJ1(x) { const ax = Math.abs(x);
  if (ax < 8) { const y = x * x; return x * (72362614232 + y * (-7895059235 + y * (242396853.1 + y * (-2972611.439 + y * (15704.4826 + y * -30.16036606))))) / (144725228442 + y * (2300535178 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))))); }
  const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
  const v = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * (1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)))) - z * Math.sin(xx) * (0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)))));
  return x < 0 ? -v : v; }
// Zeros of J1 (after the one at 0): Newton steps from (n + 1/4)·π.
const J1ZEROS = Array.from({ length: 60 }, (_, i) => { let a = (i + 1.25) * Math.PI; for (let k = 0; k < 6; k++) a -= besselJ1(a) / (besselJ0(a) - besselJ1(a) / a); return a; });
// Steady temperature in a disc of radius R and thickness H, heated on top by a Gaussian flux q0·exp(−2r²/w²),
// insulated on its side and held at 0 on its bottom face. Returns T(r, z) in kelvin above the bottom face.
function discSolver(q0, w, R, H, k) {
  const M = 400, q = r => q0 * Math.exp(-2 * r * r / (w * w));
  const integ = f => { let s = 0; for (let i = 0; i < M; i++) { const r = (i + .5) * R / M; s += f(r) * r; } return s * R / M; };
  const c0 = 2 / (R * R) * integ(q);
  const cs = J1ZEROS.map(a => { const j = besselJ0(a); return 2 / (R * R * j * j) * integ(r => q(r) * besselJ0(a * r / R)); });
  return (r, z) => { let T = c0 * (H - z) / k;
    J1ZEROS.forEach((a, n) => { const e = Math.exp(-2 * a * H / R), u = a * (H - z) / R, t = a * z / R;
      // sinh(u)/cosh(aH/R), written with decaying exponentials only
      T += cs[n] * R / (k * a) * besselJ0(a * r / R) * (Math.exp(-t) - Math.exp(-t - 2 * u)) / (1 + e); });
    return T; };
}

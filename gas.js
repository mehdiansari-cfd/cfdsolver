/* ============================================================
   gas.js  —  gas models in SI units (math layer, no DOM)
   ------------------------------------------------------------
   Two thermodynamic models, both ideal gases (p = rho R T):

     calorically perfect  cp constant          → gamma constant
     thermally perfect    cp = cp(T) (NASA-7)  → gamma = gamma(T)

   Species data: NASA 7-coefficient polynomials, 200–6000 K, from
     B.J. McBride, S. Gordon, M.A. Reno, "Coefficients for Calculating
     Thermodynamic and Transport Properties of Individual Species",
     NASA TM-4513 (1993), as distributed in Cantera's nasa_gas.yaml.
   Outside 200–6000 K, cp is held at its end value (constant-cp
   extension) so h, phi and their inverses stay smooth and monotone.

   A gas model is a plain object, per unit mass, in SI:
     { name, label, kind: "thermally-perfect" | "calorically-perfect",
       R      [J/(kg K)]  specific gas constant
       cp(T)  [J/(kg K)]
       h(T)   [J/kg]      sensible enthalpy, h(Tref) = 0
       phi(T) [J/(kg K)]  entropy function ∫ cp/T dT, phi(Tref) = 0
                          (s = phi(T) - R ln(p/pref) + const) }
   Everything else (e, gamma, a, T from e, isentropes, stagnation
   states) is derived from those four by the helpers below, so they
   work for every model, including binary mixtures (mixGases).

   Public API (browser globals + CommonJS):
     GASES, gasModel(name, model), caloricallyPerfect(gas, T)
     mixGases(A, B, Y)            mass fraction Y of A
     gasE, gasGamma, gasSound, gasTfromE, gasTfromH, gasStagnation,
     gasStaticFromStagnation, gasIsentropeT, RU, TREF, PREF
   ============================================================ */

const RU   = 8.314462618;   // universal gas constant [J/(mol K)] (CODATA 2018)
const TREF = 298.15;        // thermochemical reference temperature [K]
const PREF = 101325;        // reference pressure [Pa]
const NEWTON_TOL = 1e-10;   // relative tolerance on T in the inverse solves
const NEWTON_MAX = 50;

// ---------- Species: NASA-7 coefficients + molar mass ----------
// lo: [Tmin, Tmid], hi: [Tmid, Tmax]. W in kg/mol (IUPAC standard atomic weights).
const SPECIES = {
  N2:  { W: 28.0134e-3, T: [200, 1000, 6000],
         lo: [3.53100528, -1.23660987e-04, -5.02999437e-07, 2.43530612e-09, -1.40881235e-12, -1046.97628, 2.96747468],
         hi: [2.95257626, 1.39690057e-03, -4.92631691e-07, 7.86010367e-11, -4.60755321e-15, -923.948645, 5.87189252] },
  O2:  { W: 31.9988e-3, T: [200, 1000, 6000],
         lo: [3.78245636, -2.99673415e-03, 9.847302e-06, -9.68129508e-09, 3.24372836e-12, -1063.94356, 3.65767573],
         hi: [3.66096083, 6.56365523e-04, -1.41149485e-07, 2.05797658e-11, -1.29913248e-15, -1215.97725, 3.41536184] },
  Ar:  { W: 39.948e-3, T: [200, 6000, 6000],
         lo: [2.5, 0, 0, 0, 0, -745.375, 4.37967491],
         hi: [2.5, 0, 0, 0, 0, -745.375, 4.37967491] },
  He:  { W: 4.002602e-3, T: [200, 6000, 6000],
         lo: [2.5, 0, 0, 0, 0, -745.375, 0.928724724],
         hi: [2.5, 0, 0, 0, 0, -745.375, 0.928724724] },
  H2:  { W: 2.01588e-3, T: [200, 1000, 6000],
         lo: [2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12, -917.935173, 0.683010238],
         hi: [2.93286579, 8.26607967e-04, -1.46402335e-07, 1.54100359e-11, -6.88804432e-16, -813.065597, -1.02432887] },
  H2O: { W: 18.01528e-3, T: [200, 1000, 6000],
         lo: [4.19864056, -2.0364341e-03, 6.52040211e-06, -5.48797062e-09, 1.77197817e-12, -3.02937267e+04, -0.849032208],
         hi: [2.67703787, 2.97318329e-03, -7.7376969e-07, 9.44336689e-11, -4.26900959e-15, -2.98858938e+04, 6.88255571] },
  CO2: { W: 44.0095e-3, T: [200, 1000, 6000],
         lo: [2.35677352, 8.98459677e-03, -7.12356269e-06, 2.45919022e-09, -1.43699548e-13, -4.83719697e+04, 9.90105222],
         hi: [4.63659493, 2.74131991e-03, -9.95828531e-07, 1.60373011e-10, -9.16103468e-15, -4.90249341e+04, -1.93534855] },
  CO:  { W: 28.0101e-3, T: [200, 1000, 6000],
         lo: [3.57953347, -6.1035368e-04, 1.01681433e-06, 9.07005884e-10, -9.04424499e-13, -1.4344086e+04, 3.50840928],
         hi: [3.04848583, 1.35172818e-03, -4.85794075e-07, 7.88536486e-11, -4.69807489e-15, -1.42661171e+04, 6.0170979] }
};

// NASA-7 dimensionless forms (McBride et al. Eqs. 1–3):
//   cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
//   h/RT = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
//   s/R  = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7
function _nasaCoeffs(sp, T) { return T < sp.T[1] ? sp.lo : sp.hi; }
function _cpR(a, T) { return a[0] + T * (a[1] + T * (a[2] + T * (a[3] + T * a[4]))); }
function _hR(a, T)  { return a[5] + T * (a[0] + T * (a[1] / 2 + T * (a[2] / 3 + T * (a[3] / 4 + T * a[4] / 5)))); }
function _sR(a, T)  { return a[0] * Math.log(T) + a[6] + T * (a[1] + T * (a[2] / 2 + T * (a[3] / 3 + T * a[4] / 4))); }

// Species functions in molar-dimensionless form with the constant-cp
// extension outside the fitted range.
function _speciesCpR(sp, T) {
  const Tc = Math.min(sp.T[2], Math.max(sp.T[0], T));
  return _cpR(_nasaCoeffs(sp, Tc), Tc);
}
function _speciesHR(sp, T) {          // h / R  [K]
  const Tc = Math.min(sp.T[2], Math.max(sp.T[0], T));
  const a = _nasaCoeffs(sp, Tc);
  return _hR(a, Tc) + _cpR(a, Tc) * (T - Tc);
}
function _speciesSR(sp, T) {          // s°(T) / R  [-]
  const Tc = Math.min(sp.T[2], Math.max(sp.T[0], T));
  const a = _nasaCoeffs(sp, Tc);
  return _sR(a, Tc) + _cpR(a, Tc) * Math.log(T / Tc);
}

// ---------- Thermally perfect mixture of fixed composition ----------
// X: mole fractions { species: x }. Returns a per-unit-mass gas model.
function _mixtureFromMoles(name, label, X, note) {
  const keys = Object.keys(X);
  let xs = 0, W = 0;
  for (const k of keys) xs += X[k];
  for (const k of keys) W += (X[k] / xs) * SPECIES[k].W;
  const x = keys.map(k => X[k] / xs), sp = keys.map(k => SPECIES[k]);
  const R = RU / W;
  // Mixture molar quantities are mole-fraction sums; per unit mass divide by W.
  const cpR = T => { let s = 0; for (let i = 0; i < sp.length; i++) s += x[i] * _speciesCpR(sp[i], T); return s; };
  const hR  = T => { let s = 0; for (let i = 0; i < sp.length; i++) s += x[i] * _speciesHR(sp[i], T); return s; };
  const sR  = T => { let s = 0; for (let i = 0; i < sp.length; i++) s += x[i] * _speciesSR(sp[i], T); return s; };
  const hRef = hR(TREF), sRef = sR(TREF);
  return {
    name, label, note: note || "", kind: "thermally-perfect", R, W, composition: X,
    cp:  T => R * cpR(T),
    h:   T => R * (hR(T) - hRef),            // sensible enthalpy, h(Tref) = 0
    phi: T => R * (sR(T) - sRef)             // ∫ cp/T dT from Tref
  };
}

// ---------- Catalogue ----------
// Combustion products are frozen at complete-combustion, stoichiometric
// composition (no dissociation): transparent and reproducible, but they
// over-predict flame temperature. Fine for gamma(T) and R; not for chemistry.
const GAS_DEFS = {
  air:  { label: "Air (dry)",          X: { N2: 0.78084, O2: 0.20946, Ar: 0.00934, CO2: 0.00036 } },
  N2:   { label: "Nitrogen N₂",        X: { N2: 1 } },
  O2:   { label: "Oxygen O₂",          X: { O2: 1 } },
  Ar:   { label: "Argon Ar",           X: { Ar: 1 } },
  He:   { label: "Helium He",          X: { He: 1 } },
  H2:   { label: "Hydrogen H₂",        X: { H2: 1 } },
  CO2:  { label: "Carbon dioxide CO₂", X: { CO2: 1 } },
  "CH4-air products": { label: "CH₄–air products (φ = 1, frozen)", X: { CO2: 1, H2O: 2, N2: 7.52 },
                        note: "CH4 + 2(O2 + 3.76 N2) → CO2 + 2 H2O + 7.52 N2, no dissociation" },
  "H2-air products":  { label: "H₂–air products (φ = 1, frozen)",  X: { H2O: 1, N2: 1.88 },
                        note: "H2 + ½(O2 + 3.76 N2) → H2O + 1.88 N2, no dissociation" }
};
const GASES = {};
for (const k of Object.keys(GAS_DEFS)) {
  GASES[k] = _mixtureFromMoles(k, GAS_DEFS[k].label, GAS_DEFS[k].X, GAS_DEFS[k].note);
}

// Calorically perfect version of a gas: cp frozen at temperature T (default Tref).
function caloricallyPerfect(gas, T) {
  const Tf = T || TREF, cp = gas.cp(Tf);
  return {
    name: gas.name, label: gas.label, note: gas.note, kind: "calorically-perfect",
    R: gas.R, W: gas.W, composition: gas.composition, cpFrozenAt: Tf,
    cp:  () => cp,
    h:   T => cp * (T - TREF),
    phi: T => cp * Math.log(T / TREF)
  };
}

// Look up a gas by name; model = "thermally-perfect" (default) | "calorically-perfect".
function gasModel(name, model) {
  const g = GASES[name];
  if (!g) throw new Error("Unknown gas '" + name + "'. Known: " + Object.keys(GASES).join(", "));
  return model === "calorically-perfect" ? caloricallyPerfect(g) : g;
}

// Binary mixture by mass fraction Y of gas A (1 - Y of gas B), e.g. the
// driver / driven gases mixing across a numerically smeared contact.
// Ideal-gas mixing: per-mass properties are mass-fraction weighted.
function mixGases(A, B, Y) {
  if (Y >= 1) return A;
  if (Y <= 0) return B;
  const Z = 1 - Y;
  return {
    name: A.name + "+" + B.name, kind: A.kind === B.kind ? A.kind : "thermally-perfect",
    R: Y * A.R + Z * B.R,
    cp:  T => Y * A.cp(T) + Z * B.cp(T),
    h:   T => Y * A.h(T) + Z * B.h(T),
    phi: T => Y * A.phi(T) + Z * B.phi(T)
  };
}

// ---------- Derived properties (work for every model) ----------
function gasE(gas, T)     { return gas.h(T) - gas.R * T; }                      // e = h - RT
function gasGamma(gas, T) { const cp = gas.cp(T); return cp / (cp - gas.R); }
function gasSound(gas, T) { return Math.sqrt(gasGamma(gas, T) * gas.R * T); }   // a^2 = gamma R T

// Newton inverse of a monotone function f(T) with derivative fp(T) > 0.
function _invert(f, fp, target, guess) {
  let T = guess > 0 ? guess : TREF;
  for (let it = 0; it < NEWTON_MAX; it++) {
    const dT = (f(T) - target) / fp(T);
    let Tn = T - dT;
    if (!(Tn > 0)) Tn = 0.5 * T;                 // stay on T > 0
    if (Math.abs(Tn - T) <= NEWTON_TOL * T) return Tn;
    T = Tn;
  }
  return T;
}
// T from specific internal energy: de/dT = cv = cp - R.
function gasTfromE(gas, e, guess) {
  return _invert(T => gasE(gas, T), T => gas.cp(T) - gas.R, e, guess);
}
// T from specific enthalpy: dh/dT = cp.
function gasTfromH(gas, h, guess) {
  return _invert(gas.h, gas.cp, h, guess);
}
// T on the isentrope through (T1, p1) at pressure p: phi(T) - R ln p = const.
function gasIsentropeT(gas, T1, p1, p) {
  const target = gas.phi(T1) + gas.R * Math.log(p / p1);
  return _invert(gas.phi, T => gas.cp(T) / T, target, T1 * Math.pow(p / p1, gas.R / gas.cp(T1)));
}
// Stagnation state of (T, p, u): h0 = h + u^2/2, isentropic p0.
function gasStagnation(gas, T, p, u) {
  const T0 = gasTfromH(gas, gas.h(T) + 0.5 * u * u, T);
  return { T0, p0: p * Math.exp((gas.phi(T0) - gas.phi(T)) / gas.R) };
}
// Static state at speed u from stagnation (T0, p0).
function gasStaticFromStagnation(gas, T0, p0, u) {
  const T = gasTfromH(gas, gas.h(T0) - 0.5 * u * u, T0);
  return { T, p: p0 * Math.exp((gas.phi(T) - gas.phi(T0)) / gas.R) };
}

if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    RU, TREF, PREF, SPECIES, GASES, GAS_DEFS,
    gasModel, caloricallyPerfect, mixGases,
    gasE, gasGamma, gasSound, gasTfromE, gasTfromH, gasIsentropeT,
    gasStagnation, gasStaticFromStagnation
  };
}

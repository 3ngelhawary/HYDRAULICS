// Hydraulic Calculator - Pressure + Gravity (Manning + partial flow d/D)
// No dependencies.

const $ = (id) => document.getElementById(id);
function clamp(v, a, b){ return Math.max(a, Math.min(b, v)); }

function toM3s(qValue, unit){
  const q = Number(qValue);
  if (!isFinite(q)) return NaN;

  switch(unit){
    case "Lps":  return q / 1000.0;
    case "m3ps": return q;
    case "Lpm":  return q / 1000.0 / 60.0;
    case "m3ph": return q / 3600.0;
    case "gpm":  return q * 0.003785411784 / 60.0; // US gallon/min -> m3/s
    case "cfs":  return q * 0.028316846592;        // ft3/s -> m3/s
    default:     return NaN;
  }
}

function velToMps(v, unit){
  const x = Number(v);
  if (!isFinite(x)) return NaN;
  return unit === "ftps" ? x * 0.3048 : x;
}

function diameterToM(d, unit){
  const x = Number(d);
  if (!isFinite(x)) return NaN;
  if (unit === "mm") return x / 1000.0;
  if (unit === "in") return x * 0.0254;
  return x;
}

function lenToM(L, unit){
  const x = Number(L);
  if (!isFinite(x)) return NaN;
  if (unit === "km") return x * 1000.0;
  if (unit === "ft") return x * 0.3048;
  return x; // m
}

function slopeToMmPerM(val, unit){
  const s = Number(val);
  if (!isFinite(s)) return NaN;
  if (unit === "percent") return s / 100.0;
  if (unit === "permil")  return s / 1000.0;
  return s; // m/m
}

function formatSig(x, digits=4){
  if (!isFinite(x)) return "—";
  const abs = Math.abs(x);
  if (abs === 0) return "0";
  if (abs < 0.001 || abs >= 10000) return x.toExponential(3);
  return Number(x.toPrecision(digits)).toString();
}

function areaFromD(Dm){
  return Math.PI * Dm * Dm / 4.0;
}
function velocityFromQandD(Qm3s, Dm){
  const A = areaFromD(Dm);
  if (!isFinite(A) || A <= 0) return NaN;
  return Qm3s / A;
}
function diameterFromQandV(Qm3s, Vmps){
  if (!isFinite(Qm3s) || !isFinite(Vmps) || Qm3s <= 0 || Vmps <= 0) return NaN;
  return Math.sqrt((4.0 * Qm3s) / (Math.PI * Vmps));
}

function mpsToFtps(v){ return v / 0.3048; }
function mToMm(m){ return m * 1000.0; }
function mToIn(m){ return m / 0.0254; }

// Hazen–Williams (SI) headloss in meters:
// hf = 10.67 * L * Q^1.852 / (C^1.852 * D^4.871)
function hazenWilliamsHf(Lm, Qm3s, Cm, Dm){
  const C = Number(Cm);
  if (!isFinite(Lm) || !isFinite(Qm3s) || !isFinite(C) || !isFinite(Dm)) return NaN;
  if (Lm <= 0 || Qm3s <= 0 || C <= 0 || Dm <= 0) return NaN;
  return 10.67 * Lm * Math.pow(Qm3s, 1.852) / (Math.pow(C, 1.852) * Math.pow(Dm, 4.871));
}

// Manning full-flow
function manningQ_full(Dm, n, S){
  const A = areaFromD(Dm);
  const R = Dm / 4.0;
  if (!isFinite(A) || !isFinite(R) || A <= 0 || R <= 0 || n <= 0 || S <= 0) return NaN;
  return (1.0 / n) * A * Math.pow(R, 2.0/3.0) * Math.sqrt(S);
}

// Solve D for full-flow Manning
function solveD_manning_full(Qtarget, n, S){
  let lo = 0.01;
  let hi = 5.0;
  if (!isFinite(Qtarget) || !isFinite(n) || !isFinite(S) || Qtarget <= 0 || n <= 0 || S <= 0) return NaN;

  for (let i=0; i<40; i++){
    const qHi = manningQ_full(hi, n, S);
    if (isFinite(qHi) && qHi >= Qtarget) break;
    hi *= 1.4;
    if (hi > 50) break;
  }

  for (let iter=0; iter<80; iter++){
    const mid = (lo + hi) / 2.0;
    const qMid = manningQ_full(mid, n, S);
    if (!isFinite(qMid)) return NaN;
    if (qMid < Qtarget) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2.0;
}

// Partial flow geometry for circular conduit (open channel):
// theta in [0, 2π] (wetted angle).
// A = (D^2/8) * (theta - sin(theta))
// P = (D/2) * theta
// R = A / P
// depth ratio y/D = (1 - cos(theta/2)) / 2
function manningQ_partial(Dm, n, S, theta){
  if (!isFinite(Dm) || !isFinite(n) || !isFinite(S) || !isFinite(theta)) return NaN;
  if (Dm <= 0 || n <= 0 || S <= 0) return NaN;

  const A = (Dm*Dm/8.0) * (theta - Math.sin(theta));
  const P = (Dm/2.0) * theta;

  if (!isFinite(A) || !isFinite(P) || A <= 0 || P <= 0) return NaN;

  const R = A / P;
  return (1.0 / n) * A * Math.pow(R, 2.0/3.0) * Math.sqrt(S);
}

function depthRatioFromTheta(theta){
  // y/D
  return (1.0 - Math.cos(theta/2.0)) / 2.0;
}

// Solve theta for given Q with fixed D (returns theta, y/D, Qfull)
function solveThetaForQ(Dm, Qtarget, n, S){
  if (!isFinite(Dm) || !isFinite(Qtarget) || !isFinite(n) || !isFinite(S)) return { theta: NaN, yd: NaN, qfull: NaN };
  if (Dm <= 0 || Qtarget <= 0 || n <= 0 || S <= 0) return { theta: NaN, yd: NaN, qfull: NaN };

  const qfull = manningQ_full(Dm, n, S);
  if (!isFinite(qfull) || qfull <= 0) return { theta: NaN, yd: NaN, qfull: NaN };

  // If flow exceeds full capacity -> cap at 100%
  if (Qtarget >= qfull) {
    return { theta: 2*Math.PI, yd: 1.0, qfull };
  }

  let lo = 1e-6;          // near zero
  let hi = 2*Math.PI;     // full

  // bisection
  for (let iter=0; iter<90; iter++){
    const mid = (lo + hi) / 2.0;
    const qMid = manningQ_partial(Dm, n, S, mid);
    if (!isFinite(qMid)) return { theta: NaN, yd: NaN, qfull: NaN };
    if (qMid < Qtarget) lo = mid; else hi = mid;
  }

  const theta = (lo + hi) / 2.0;
  const yd = depthRatioFromTheta(theta);
  return { theta, yd, qfull };
}

// UI helpers
function setLive(mode, bigNumber, bigUnit, fillPct, hint){
  $("modePill").textContent = mode;
  $("bigNumber").textContent = bigNumber;
  $("bigUnit").textContent = bigUnit;
  $("meterFill").style.width = `${clamp(fillPct, 0, 100)}%`;
  $("meterHint").textContent = hint;
}

function gaugeForVelocity(vMps){
  const pct = clamp((vMps / 4.0) * 100.0, 0, 100);
  $("gaugeNeedle").style.left = `${pct}%`;

  if (!isFinite(vMps)) { $("gaugePill").textContent = "—"; return; }
  if (vMps < 1.0) $("gaugePill").textContent = "Low";
  else if (vMps <= 2.5) $("gaugePill").textContent = "OK";
  else if (vMps <= 3.5) $("gaugePill").textContent = "High";
  else $("gaugePill").textContent = "Very High";
}

function updatePressureModeUI(){
  const mode = $("pressureMode").value;
  const showTarget = (mode === "D_from_V");
  $("targetVWrap").style.display = showTarget ? "" : "none";
  $("diamWrap").style.display = showTarget ? "none" : "";
}

function updateGravityModeUI(){
  const mode = $("gravityMode").value;
  // If finding diameter, we don't need user diameter input
  $("gDiamWrap").style.display = (mode === "DD_from_D") ? "" : "none";
}

function calcPressure(){
  const count = Math.max(1, Math.floor(Number($("pCount").value) || 1));
  const Q = toM3s($("qVal").value, $("qUnit").value);
  const mode = $("pressureMode").value;

  let Dm = NaN;
  let Vmps = NaN;

  if (mode === "D_from_V"){
    const Vt = velToMps($("vTarget").value, $("vUnit").value);
    Dm = diameterFromQandV(Q, Vt);
    Vmps = Vt;
  } else {
    Dm = diameterToM($("dVal").value, $("dUnit").value);
    Vmps = velocityFromQandD(Q, Dm);
  }

  const Lm = lenToM($("pLen").value, $("pLenUnit").value);
  const C = Number($("hwC").value);
  const hf = hazenWilliamsHf(Lm, Q, C, Dm);
  const hfTot = isFinite(hf) ? hf * count : NaN;

  // Results
  $("resD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("resDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 5)} m` : "—";

  $("resV").textContent = isFinite(Vmps) ? `${formatSig(Vmps, 5)} m/s` : "—";
  $("resVSub").textContent = isFinite(Vmps) ? `${formatSig(mpsToFtps(Vmps), 5)} ft/s` : "—";

  $("resHf").textContent = isFinite(hf) ? `${formatSig(hf, 6)} m` : "—";
  $("resHfSub").textContent = (isFinite(hf) && isFinite(Lm) && Lm>0)
    ? `${formatSig(hf / (Lm/1000.0), 6)} m/km`
    : "—";

  $("resHfTot").textContent = isFinite(hfTot) ? `${formatSig(hfTot, 6)} m` : "—";
  $("resHfTotSub").textContent = isFinite(hfTot) ? `${count} pipe(s)` : "—";

  // Quick check
  $("miniQ").textContent = isFinite(Q) ? formatSig(Q, 6) : "—";
  $("miniD").textContent = isFinite(Dm) ? formatSig(Dm, 6) : "—";
  $("miniV").textContent = isFinite(Vmps) ? formatSig(Vmps, 6) : "—";
  $("miniHf").textContent = isFinite(hf) ? formatSig(hf, 6) : "—";

  // Live: show total headloss (more useful when multiple pipes)
  const big = isFinite(hfTot) ? formatSig(hfTot, 6) : "—";
  setLive("Pressure", big, "m headloss", isFinite(hfTot) ? clamp((hfTot/50)*100,0,100) : 0, "Total Hazen–Williams headloss.");

  gaugeForVelocity(Vmps);
}

function resetPressure(){
  $("pCount").value = "1";
  $("qVal").value = "";
  $("qUnit").value = "Lps";
  $("pressureMode").value = "D_from_V";
  $("vTarget").value = "";
  $("vUnit").value = "mps";
  $("dVal").value = "";
  $("dUnit").value = "mm";
  $("pLen").value = "";
  $("pLenUnit").value = "m";
  $("hwC").value = "130";

  updatePressureModeUI();

  $("resD").textContent = "—";
  $("resDSub").textContent = "—";
  $("resV").textContent = "—";
  $("resVSub").textContent = "—";
  $("resHf").textContent = "—";
  $("resHfSub").textContent = "—";
  $("resHfTot").textContent = "—";
  $("resHfTotSub").textContent = "—";

  $("miniQ").textContent = "—";
  $("miniD").textContent = "—";
  $("miniV").textContent = "—";
  $("miniHf").textContent = "—";

  setLive("Pressure", "—", "—", 0, "Enter inputs to see results.");
  gaugeForVelocity(NaN);
}

function calcGravity(){
  const Q = toM3s($("gQVal").value, $("gQUnit").value);
  const n = Number($("nVal").value);
  const S = slopeToMmPerM($("sVal").value, $("sUnit").value);
  const mode = $("gravityMode").value;

  let Dm = NaN;
  let yd = NaN;
  let qfull = NaN;

  if (mode === "D_full_from_Q"){
    Dm = solveD_manning_full(Q, n, S);
    // If we sized for full flow, d/D is ~100%
    yd = isFinite(Dm) ? 1.0 : NaN;
    qfull = isFinite(Dm) ? manningQ_full(Dm, n, S) : NaN;
  } else {
    Dm = diameterToM($("gDVal").value, $("gDUnit").value);
    const solved = solveThetaForQ(Dm, Q, n, S);
    yd = solved.yd;
    qfull = solved.qfull;
  }

  const Vmps = velocityFromQandD(Q, Dm);

  $("gResD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("gResDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 5)} m` : "—";

  $("gResV").textContent = isFinite(Vmps) ? `${formatSig(Vmps, 5)} m/s` : "—";
  $("gResVSub").textContent = isFinite(Vmps) ? `${formatSig(mpsToFtps(Vmps), 5)} ft/s` : "—";

  $("gResDD").textContent = isFinite(yd) ? `${formatSig(yd*100.0, 5)} %` : "—";
  $("gResDDSub").textContent = isFinite(qfull) && isFinite(Q)
    ? `Q / Qfull = ${formatSig((Q/qfull)*100.0, 5)} %`
    : "—";

  $("gResQfull").textContent = isFinite(qfull) ? `${formatSig(qfull, 6)} m³/s` : "—";

  $("gMiniQ").textContent = isFinite(Q) ? formatSig(Q, 6) : "—";
  $("gMiniS").textContent = isFinite(S) ? formatSig(S, 6) : "—";
  $("gMiniN").textContent = isFinite(n) ? formatSig(n, 6) : "—";
  $("gMiniDD").textContent = isFinite(yd) ? formatSig(yd, 6) : "—";

  // Live: show d/D %
  const big = isFinite(yd) ? formatSig(yd*100.0, 5) : "—";
  setLive("Gravity", big, "% d/D", isFinite(yd) ? clamp(yd*100,0,100) : 0, "Depth ratio (Manning partial flow).");
}

function resetGravity(){
  $("gCount").value = "1";
  $("gQVal").value = "";
  $("gQUnit").value = "Lps";
  $("nVal").value = "0.013";
  $("sVal").value = "";
  $("sUnit").value = "mperm";
  $("gravityMode").value = "DD_from_D";
  $("gDVal").value = "";
  $("gDUnit").value = "mm";

  updateGravityModeUI();

  $("gResD").textContent = "—";
  $("gResDSub").textContent = "—";
  $("gResV").textContent = "—";
  $("gResVSub").textContent = "—";
  $("gResDD").textContent = "—";
  $("gResDDSub").textContent = "—";
  $("gResQfull").textContent = "—";

  $("gMiniQ").textContent = "—";
  $("gMiniS").textContent = "—";
  $("gMiniN").textContent = "—";
  $("gMiniDD").textContent = "—";

  setLive("Gravity", "—", "—", 0, "Enter inputs to see results.");
}

function setTab(tabName){
  document.querySelectorAll(".tab").forEach(btn=>{
    const active = btn.dataset.tab === tabName;
    btn.classList.toggle("is-active", active);
    btn.setAttribute("aria-selected", active ? "true" : "false");
  });
  document.querySelectorAll(".tabpane").forEach(p=>{
    p.classList.toggle("is-active", p.dataset.pane === tabName);
  });

  // Keep whatever is currently shown; next calc will update mode label correctly
  $("modePill").textContent = tabName === "pressure" ? "Pressure" : "Gravity";
}

// Wire up
document.addEventListener("DOMContentLoaded", ()=>{
  // Tabs
  document.querySelectorAll(".tab").forEach(btn=>{
    btn.addEventListener("click", ()=> setTab(btn.dataset.tab));
  });

  // Pressure
  $("pressureMode").addEventListener("change", updatePressureModeUI);
  updatePressureModeUI();

  $("btnCalcPressure").addEventListener("click", calcPressure);
  $("btnResetPressure").addEventListener("click", resetPressure);

  ["pCount","qVal","qUnit","vTarget","vUnit","dVal","dUnit","pLen","pLenUnit","hwC","pressureMode"].forEach(id=>{
    $(id).addEventListener("input", ()=> {
      if ($("qVal").value || $("vTarget").value || $("dVal").value || $("pLen").value) calcPressure();
    });
    $(id).addEventListener("change", ()=> {
      if ($("qVal").value || $("vTarget").value || $("dVal").value || $("pLen").value) calcPressure();
    });
  });

  // Gravity
  $("gravityMode").addEventListener("change", updateGravityModeUI);
  updateGravityModeUI();

  $("btnCalcGravity").addEventListener("click", calcGravity);
  $("btnResetGravity").addEventListener("click", resetGravity);

  ["gCount","gQVal","gQUnit","nVal","sVal","sUnit","gravityMode","gDVal","gDUnit"].forEach(id=>{
    $(id).addEventListener("input", ()=> {
      if ($("gQVal").value || $("sVal").value || $("gDVal").value) calcGravity();
    });
    $(id).addEventListener("change", ()=> {
      if ($("gQVal").value || $("sVal").value || $("gDVal").value) calcGravity();
    });
  });

  // Defaults
  resetPressure();
});

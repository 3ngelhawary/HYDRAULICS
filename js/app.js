// Hydraulic Calculator - Pressure + Gravity (Manning + partial flow d/D)
// Total flow input; per-pipe flow = total / number of parallel pipes.

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
    case "gpm":  return q * 0.003785411784 / 60.0;
    case "cfs":  return q * 0.028316846592;
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
  return x;
}

function slopeToMmPerM(val, unit){
  const s = Number(val);
  if (!isFinite(s)) return NaN;
  if (unit === "percent") return s / 100.0;
  if (unit === "permil")  return s / 1000.0;
  return s;
}

function formatSig(x, digits=4){
  if (!isFinite(x)) return "—";
  const abs = Math.abs(x);
  if (abs === 0) return "0";
  if (abs < 0.001 || abs >= 10000) return x.toExponential(3);
  return Number(x.toPrecision(digits)).toString();
}

function areaFullFromD(Dm){ return Math.PI * Dm * Dm / 4.0; }

function diameterFromQandV(Qm3s, Vmps){
  if (!isFinite(Qm3s) || !isFinite(Vmps) || Qm3s <= 0 || Vmps <= 0) return NaN;
  return Math.sqrt((4.0 * Qm3s) / (Math.PI * Vmps));
}

function velocityFull(Qm3s, Dm){
  const A = areaFullFromD(Dm);
  if (!isFinite(A) || A <= 0) return NaN;
  return Qm3s / A;
}

function mpsToFtps(v){ return v / 0.3048; }
function mToMm(m){ return m * 1000.0; }
function mToIn(m){ return m / 0.0254; }

// Hazen–Williams (SI) headloss in meters (Q per pipe):
// hf = 10.67 * L * Q^1.852 / (C^1.852 * D^4.871)
function hazenWilliamsHf(Lm, Qm3s, C, Dm){
  const c = Number(C);
  if (!isFinite(Lm) || !isFinite(Qm3s) || !isFinite(c) || !isFinite(Dm)) return NaN;
  if (Lm <= 0 || Qm3s <= 0 || c <= 0 || Dm <= 0) return NaN;
  return 10.67 * Lm * Math.pow(Qm3s, 1.852) / (Math.pow(c, 1.852) * Math.pow(Dm, 4.871));
}

// Manning full-flow circular
function manningQ_full(Dm, n, S){
  const A = areaFullFromD(Dm);
  const R = Dm / 4.0;
  if (!isFinite(A) || !isFinite(R) || A <= 0 || R <= 0 || n <= 0 || S <= 0) return NaN;
  return (1.0 / n) * A * Math.pow(R, 2.0/3.0) * Math.sqrt(S);
}

// Partial flow (open channel) circular geometry using wetted angle theta [0..2π]:
// A = (D^2/8) * (theta - sin(theta))
// P = (D/2) * theta
// R = A / P
// y/D = (1 - cos(theta/2)) / 2
function areaPartial(Dm, theta){
  return (Dm*Dm/8.0) * (theta - Math.sin(theta));
}
function manningQ_partial(Dm, n, S, theta){
  if (Dm <= 0 || n <= 0 || S <= 0) return NaN;
  const A = areaPartial(Dm, theta);
  const P = (Dm/2.0) * theta;
  if (!isFinite(A) || !isFinite(P) || A <= 0 || P <= 0) return NaN;
  const R = A / P;
  return (1.0 / n) * A * Math.pow(R, 2.0/3.0) * Math.sqrt(S);
}
function depthRatioFromTheta(theta){
  return (1.0 - Math.cos(theta/2.0)) / 2.0;
}

// Solve theta for Q (given D,n,S). Returns theta, y/D, A_wetted, Qfull
function solveThetaForQ(Dm, Qtarget, n, S){
  const qfull = manningQ_full(Dm, n, S);
  if (!isFinite(qfull) || qfull <= 0) return { theta: NaN, yd: NaN, A: NaN, qfull: NaN };

  if (Qtarget >= qfull) {
    const theta = 2*Math.PI;
    const A = areaPartial(Dm, theta);
    return { theta, yd: 1.0, A, qfull };
  }
  if (Qtarget <= 0) return { theta: NaN, yd: NaN, A: NaN, qfull };

  let lo = 1e-6;
  let hi = 2*Math.PI;

  for (let iter=0; iter<120; iter++){
    const mid = (lo + hi) / 2.0;
    const qMid = manningQ_partial(Dm, n, S, mid);
    if (!isFinite(qMid)) return { theta: NaN, yd: NaN, A: NaN, qfull: NaN };
    if (qMid < Qtarget) lo = mid; else hi = mid;
  }

  const theta = (lo + hi) / 2.0;
  const yd = depthRatioFromTheta(theta);
  const A = areaPartial(Dm, theta);
  return { theta, yd, A, qfull };
}

// Given desired depth ratio yd (0..1), compute theta and Q for a given D,n,S
function thetaFromDepthRatio(yd){
  const y = Number(yd);
  if (!isFinite(y) || y <= 0) return NaN;
  if (y >= 1) return 2*Math.PI;
  // y/D = (1 - cos(theta/2))/2 -> cos(theta/2) = 1 - 2y
  const c = 1.0 - 2.0*y;
  const half = Math.acos(clamp(c, -1, 1));
  return 2.0 * half; // theta
}

// Solve diameter for given Q at a specified d/D (yd), using bisection on D
function solveD_manning_atDepth(Qtarget, n, S, yd){
  if (!isFinite(Qtarget) || Qtarget <= 0 || !isFinite(n) || n <= 0 || !isFinite(S) || S <= 0) return NaN;
  if (!isFinite(yd) || yd <= 0 || yd > 1) return NaN;

  const theta = thetaFromDepthRatio(yd);
  if (!isFinite(theta) || theta <= 0) return NaN;

  let lo = 0.05;   // 50 mm
  let hi = 5.0;    // 5 m

  // increase hi until capacity at that depth exceeds Qtarget
  for (let i=0; i<50; i++){
    const qHi = manningQ_partial(hi, n, S, theta);
    if (isFinite(qHi) && qHi >= Qtarget) break;
    hi *= 1.35;
    if (hi > 50) break;
  }

  for (let iter=0; iter<110; iter++){
    const mid = (lo + hi) / 2.0;
    const qMid = manningQ_partial(mid, n, S, theta);
    if (!isFinite(qMid)) return NaN;
    if (qMid < Qtarget) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2.0;
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
  $("gDiamWrap").style.display = (mode === "DD_from_D") ? "" : "none";
  $("gDDWrap").style.display  = (mode === "D_from_Q_at_DD") ? "" : "none";
}

function calcPressure(){
  const pipes = Math.max(1, Math.floor(Number($("pCount").value) || 1));

  const Qtotal = toM3s($("qVal").value, $("qUnit").value);
  const Q = isFinite(Qtotal) ? (Qtotal / pipes) : NaN; // per pipe

  const mode = $("pressureMode").value;

  let Dm = NaN;
  let Vmps = NaN;

  if (mode === "D_from_V"){
    const Vt = velToMps($("vTarget").value, $("vUnit").value);
    Dm = diameterFromQandV(Q, Vt);
    Vmps = Vt;
  } else {
    Dm = diameterToM($("dVal").value, $("dUnit").value);
    Vmps = velocityFull(Q, Dm);
  }

  const Lm = lenToM($("pLen").value, $("pLenUnit").value);
  const C = Number($("hwC").value);

  const hf = hazenWilliamsHf(Lm, Q, C, Dm); // each pipe
  const hfSys = hf;                          // parallel system headloss = each pipe headloss

  $("resD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("resDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 6)} m` : "—";

  $("resV").textContent = isFinite(Vmps) ? `${formatSig(Vmps, 6)} m/s` : "—";
  $("resVSub").textContent = isFinite(Vmps) ? `${formatSig(mpsToFtps(Vmps), 6)} ft/s` : "—";

  $("resHf").textContent = isFinite(hf) ? `${formatSig(hf, 6)} m` : "—";
  $("resHfSub").textContent = (isFinite(hf) && isFinite(Lm) && Lm>0)
    ? `${formatSig(hf / (Lm/1000.0), 6)} m/km`
    : "—";

  $("resHfSys").textContent = isFinite(hfSys) ? `${formatSig(hfSys, 6)} m` : "—";
  $("resHfSysSub").textContent = `Parallel pipes: ${pipes} pipe(s)`;

  $("miniQTot").textContent = isFinite(Qtotal) ? formatSig(Qtotal, 7) : "—";
  $("miniQ").textContent = isFinite(Q) ? formatSig(Q, 7) : "—";
  $("miniD").textContent = isFinite(Dm) ? formatSig(Dm, 7) : "—";
  $("miniV").textContent = isFinite(Vmps) ? formatSig(Vmps, 7) : "—";
  $("miniHf").textContent = isFinite(hf) ? formatSig(hf, 7) : "—";

  const big = isFinite(hfSys) ? formatSig(hfSys, 6) : "—";
  setLive("Pressure", big, "m headloss", isFinite(hfSys) ? clamp((hfSys/50)*100,0,100) : 0, "Hazen–Williams headloss (parallel system).");

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
  $("pLen").value = "1000";
  $("pLenUnit").value = "m";
  $("hwC").value = "130";

  updatePressureModeUI();

  $("resD").textContent = "—";
  $("resDSub").textContent = "—";
  $("resV").textContent = "—";
  $("resVSub").textContent = "—";
  $("resHf").textContent = "—";
  $("resHfSub").textContent = "—";
  $("resHfSys").textContent = "—";
  $("resHfSysSub").textContent = "Parallel pipes → same headloss";

  $("miniQTot").textContent = "—";
  $("miniQ").textContent = "—";
  $("miniD").textContent = "—";
  $("miniV").textContent = "—";
  $("miniHf").textContent = "—";

  setLive("Pressure", "—", "—", 0, "Enter inputs to see results.");
  gaugeForVelocity(NaN);
}

function calcGravity(){
  const pipes = Math.max(1, Math.floor(Number($("gCount").value) || 1));

  const Qtotal = toM3s($("gQVal").value, $("gQUnit").value);
  const Q = isFinite(Qtotal) ? (Qtotal / pipes) : NaN; // per pipe

  const n = Number($("nVal").value);
  const S = slopeToMmPerM($("sVal").value, $("sUnit").value);
  const mode = $("gravityMode").value;

  let Dm = NaN;
  let yd = NaN;
  let qfull = NaN;
  let V = NaN;
  let wettedA = NaN;

  if (mode === "D_from_Q_at_DD"){
    const ddPct = Number($("gDDVal").value);
    yd = isFinite(ddPct) ? (ddPct / 100.0) : NaN;

    Dm = solveD_manning_atDepth(Q, n, S, yd);

    const theta = isFinite(yd) ? thetaFromDepthRatio(yd) : NaN;
    qfull = isFinite(Dm) ? manningQ_full(Dm, n, S) : NaN;

    if (isFinite(Dm) && isFinite(theta)){
      wettedA = areaPartial(Dm, theta);
      V = (isFinite(Q) && isFinite(wettedA) && wettedA > 0) ? (Q / wettedA) : NaN;
    }
  } else {
    Dm = diameterToM($("gDVal").value, $("gDUnit").value);
    const solved = solveThetaForQ(Dm, Q, n, S);
    yd = solved.yd;
    qfull = solved.qfull;
    wettedA = solved.A;
    V = (isFinite(Q) && isFinite(wettedA) && wettedA > 0) ? (Q / wettedA) : NaN;
  }

  $("gResD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("gResDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 6)} m` : "—";

  $("gResV").textContent = isFinite(V) ? `${formatSig(V, 6)} m/s` : "—";
  $("gResVSub").textContent = isFinite(V) ? `${formatSig(mpsToFtps(V), 6)} ft/s` : "—";

  $("gResDD").textContent = isFinite(yd) ? `${formatSig(yd*100.0, 6)} %` : "—";
  $("gResDDSub").textContent = (isFinite(Q) && isFinite(qfull) && qfull > 0)
    ? `Per-pipe Q/Qfull = ${formatSig((Q/qfull)*100.0, 6)} %  |  Pipes = ${pipes}`
    : "—";

  $("gResQfull").textContent = isFinite(qfull) ? `${formatSig(qfull, 7)} m³/s` : "—";

  $("gMiniQTot").textContent = isFinite(Qtotal) ? formatSig(Qtotal, 7) : "—";
  $("gMiniQ").textContent = isFinite(Q) ? formatSig(Q, 7) : "—";
  $("gMiniS").textContent = isFinite(S) ? formatSig(S, 7) : "—";
  $("gMiniN").textContent = isFinite(n) ? formatSig(n, 7) : "—";
  $("gMiniDD").textContent = isFinite(yd) ? formatSig(yd, 7) : "—";

  const big = isFinite(yd) ? formatSig(yd*100.0, 6) : "—";
  setLive("Gravity", big, "% d/D", isFinite(yd) ? clamp(yd*100,0,100) : 0, "Manning partial-flow depth ratio (per pipe).");
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
  $("gDDVal").value = "80";

  updateGravityModeUI();

  $("gResD").textContent = "—";
  $("gResDSub").textContent = "—";
  $("gResV").textContent = "—";
  $("gResVSub").textContent = "—";
  $("gResDD").textContent = "—";
  $("gResDDSub").textContent = "—";
  $("gResQfull").textContent = "—";

  $("gMiniQTot").textContent = "—";
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
  $("modePill").textContent = tabName === "pressure" ? "Pressure" : "Gravity";
}

document.addEventListener("DOMContentLoaded", ()=>{
  document.querySelectorAll(".tab").forEach(btn=>{
    btn.addEventListener("click", ()=> setTab(btn.dataset.tab));
  });

  // Pressure
  $("pressureMode").addEventListener("change", ()=> { updatePressureModeUI(); calcPressure(); });
  updatePressureModeUI();

  $("btnCalcPressure").addEventListener("click", calcPressure);
  $("btnResetPressure").addEventListener("click", resetPressure);

  ["pCount","qVal","qUnit","vTarget","vUnit","dVal","dUnit","pLen","pLenUnit","hwC","pressureMode"].forEach(id=>{
    $(id).addEventListener("input", ()=> calcPressure());
    $(id).addEventListener("change", ()=> calcPressure());
  });

  // Gravity
  $("gravityMode").addEventListener("change", ()=> { updateGravityModeUI(); calcGravity(); });
  updateGravityModeUI();

  $("btnCalcGravity").addEventListener("click", calcGravity);
  $("btnResetGravity").addEventListener("click", resetGravity);

  ["gCount","gQVal","gQUnit","nVal","sVal","sUnit","gravityMode","gDVal","gDUnit","gDDVal"].forEach(id=>{
    $(id).addEventListener("input", ()=> calcGravity());
    $(id).addEventListener("change", ()=> calcGravity());
  });

  resetPressure();
  resetGravity();
});

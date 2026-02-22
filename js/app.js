// Hydraulic Calculator - Pressure + Gravity (Manning Full Flow)
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
  if (unit === "ftps") return x * 0.3048;
  return x; // mps
}

function diameterToM(d, unit){
  const x = Number(d);
  if (!isFinite(x)) return NaN;
  if (unit === "mm") return x / 1000.0;
  if (unit === "in") return x * 0.0254;
  return x; // m
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

function slopeToMmPerM(val, unit){
  const s = Number(val);
  if (!isFinite(s)) return NaN;
  if (unit === "percent") return s / 100.0;
  if (unit === "permil")  return s / 1000.0;
  return s; // m/m
}

// Manning full-flow circular pipe:
// Q = (1/n) * A * R^(2/3) * S^(1/2)
// A = pi D^2 / 4, R = D/4 (for full circular)
// We'll solve for D given Q, n, S using bisection.
function manningQ_full(Dm, n, S){
  const A = areaFromD(Dm);
  const R = Dm / 4.0;
  if (!isFinite(A) || !isFinite(R) || A <= 0 || R <= 0 || n <= 0 || S <= 0) return NaN;
  return (1.0 / n) * A * Math.pow(R, 2.0/3.0) * Math.sqrt(S);
}

function solveD_manning_full(Qtarget, n, S){
  // bounds in meters
  let lo = 0.01;   // 10 mm
  let hi = 5.0;    // 5 m
  if (!isFinite(Qtarget) || !isFinite(n) || !isFinite(S) || Qtarget <= 0 || n <= 0 || S <= 0) return NaN;

  // ensure hi is enough
  for (let i=0; i<40; i++){
    const qHi = manningQ_full(hi, n, S);
    if (isFinite(qHi) && qHi >= Qtarget) break;
    hi *= 1.4;
    if (hi > 50) break;
  }

  // bisection
  for (let iter=0; iter<80; iter++){
    const mid = (lo + hi) / 2.0;
    const qMid = manningQ_full(mid, n, S);
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
  // 0..4 m/s mapped to 0..100%
  const pct = clamp((vMps / 4.0) * 100.0, 0, 100);
  $("gaugeNeedle").style.left = `${pct}%`;

  if (!isFinite(vMps)) {
    $("gaugePill").textContent = "—";
    return;
  }
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

function calcPressure(){
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

  const A = areaFromD(Dm);

  // Results
  $("resD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("resDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 5)} m` : "—";

  $("resV").textContent = isFinite(Vmps) ? `${formatSig(Vmps, 5)} m/s` : "—";
  $("resVSub").textContent = isFinite(Vmps) ? `${formatSig(mpsToFtps(Vmps), 5)} ft/s` : "—";

  $("resA").textContent = isFinite(A) ? `${formatSig(A, 5)} m²` : "—";

  $("miniQ").textContent = isFinite(Q) ? formatSig(Q, 6) : "—";
  $("miniD").textContent = isFinite(Dm) ? formatSig(Dm, 6) : "—";
  $("miniV").textContent = isFinite(Vmps) ? formatSig(Vmps, 6) : "—";

  // Live meter: show Diameter if sizing, else show Velocity
  if (mode === "D_from_V"){
    const big = isFinite(Dm) ? formatSig(mToMm(Dm), 5) : "—";
    setLive("Pressure", big, "mm", isFinite(Dm) ? clamp((mToMm(Dm)/1200)*100,0,100) : 0, "Diameter from target velocity.");
  } else {
    const big = isFinite(Vmps) ? formatSig(Vmps, 5) : "—";
    setLive("Pressure", big, "m/s", isFinite(Vmps) ? clamp((Vmps/4)*100,0,100) : 0, "Velocity from selected diameter.");
  }

  gaugeForVelocity(Vmps);
}

function resetPressure(){
  $("qVal").value = "";
  $("qUnit").value = "Lps";
  $("pressureMode").value = "D_from_V";
  $("vTarget").value = "";
  $("vUnit").value = "mps";
  $("dVal").value = "";
  $("dUnit").value = "mm";
  updatePressureModeUI();

  $("resD").textContent = "—";
  $("resDSub").textContent = "—";
  $("resV").textContent = "—";
  $("resVSub").textContent = "—";
  $("resA").textContent = "—";

  $("miniQ").textContent = "—";
  $("miniD").textContent = "—";
  $("miniV").textContent = "—";

  setLive("Pressure", "—", "—", 0, "Enter inputs to see results.");
  gaugeForVelocity(NaN);
}

function calcGravity(){
  const Q = toM3s($("gQVal").value, $("gQUnit").value);
  const n = Number($("nVal").value);
  const S = slopeToMmPerM($("sVal").value, $("sUnit").value);

  const Dm = solveD_manning_full(Q, n, S);
  const Qcheck = manningQ_full(Dm, n, S);
  const Vmps = velocityFromQandD(Q, Dm);

  $("gResD").textContent = isFinite(Dm) ? `${formatSig(mToMm(Dm), 5)} mm` : "—";
  $("gResDSub").textContent = isFinite(Dm) ? `${formatSig(mToIn(Dm), 5)} in   |   ${formatSig(Dm, 5)} m` : "—";

  $("gResV").textContent = isFinite(Vmps) ? `${formatSig(Vmps, 5)} m/s` : "—";
  $("gResVSub").textContent = isFinite(Vmps) ? `${formatSig(mpsToFtps(Vmps), 5)} ft/s` : "—";

  $("gResQ").textContent = isFinite(Qcheck) ? `${formatSig(Qcheck, 6)} m³/s` : "—";

  $("gMiniQ").textContent = isFinite(Q) ? formatSig(Q, 6) : "—";
  $("gMiniS").textContent = isFinite(S) ? formatSig(S, 6) : "—";
  $("gMiniN").textContent = isFinite(n) ? formatSig(n, 6) : "—";

  const big = isFinite(Dm) ? formatSig(mToMm(Dm), 5) : "—";
  setLive("Gravity", big, "mm", isFinite(Dm) ? clamp((mToMm(Dm)/1800)*100,0,100) : 0, "Manning (full flow) sizing.");
}

function resetGravity(){
  $("gQVal").value = "";
  $("gQUnit").value = "Lps";
  $("nVal").value = "0.013";
  $("sVal").value = "";
  $("sUnit").value = "mperm";

  $("gResD").textContent = "—";
  $("gResDSub").textContent = "—";
  $("gResV").textContent = "—";
  $("gResVSub").textContent = "—";
  $("gResQ").textContent = "—";

  $("gMiniQ").textContent = "—";
  $("gMiniS").textContent = "—";
  $("gMiniN").textContent = "—";

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

  if (tabName === "pressure") setLive("Pressure", $("bigNumber").textContent, $("bigUnit").textContent, parseFloat($("meterFill").style.width)||0, $("meterHint").textContent);
  else setLive("Gravity", $("bigNumber").textContent, $("bigUnit").textContent, parseFloat($("meterFill").style.width)||0, $("meterHint").textContent);
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

  // Recalc on input changes (smooth live feel)
  ["qVal","qUnit","vTarget","vUnit","dVal","dUnit","pressureMode"].forEach(id=>{
    $(id).addEventListener("input", ()=> {
      // Only auto-calc if something meaningful is typed
      if (($("qVal").value || $("vTarget").value || $("dVal").value)) calcPressure();
    });
    $(id).addEventListener("change", ()=> {
      if (($("qVal").value || $("vTarget").value || $("dVal").value)) calcPressure();
    });
  });

  // Gravity
  $("btnCalcGravity").addEventListener("click", calcGravity);
  $("btnResetGravity").addEventListener("click", resetGravity);

  ["gQVal","gQUnit","nVal","sVal","sUnit"].forEach(id=>{
    $(id).addEventListener("input", ()=> {
      if ($("gQVal").value || $("sVal").value) calcGravity();
    });
    $(id).addEventListener("change", ()=> {
      if ($("gQVal").value || $("sVal").value) calcGravity();
    });
  });

  // Defaults
  resetPressure();
});

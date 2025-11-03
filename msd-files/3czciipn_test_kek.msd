// =============================================================
// MultiSystem model for TADF emitters (no CBP) — 3CzClIPN default
// PP (polaron pair)  →  X (exciton: S1, T0, Tp, Tm)  →  SINK (RAD/NR)
// Units: fields in Tesla; ZFS-like principal values in mT (via prefactor=1e-3); rates in ns^-1
// Sweep: 401 steps from -0.2005 T in +0.001 T increments, applied to both PP and X Zeeman fields
// Time-evolution: MultiStaticSS-TimeEvolution (finite window; sinks accumulate)
// =============================================================

// ---------------------- PP: polaron pair (F-pair initial mixture) ----------------------
SpinSystem PP
{
  // electrons/holes with identical g (Δg ≈ 0 in PP)
  Spin e  { type = electron; spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 2.0023"); }
  Spin h  { type = electron; spin = 1/2; tensor = isotropic("0 0 0; 0 0 0; 0 0 2.0023"); }

  // two effective nuclei giving ~6 mT LF width (rms) via isotropic HF 3.8
  Spin n1 { type = nucleus;  spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 6.0"); }   // mT
  Spin n2 { type = nucleus;  spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 6.0"); }   // mT

  // interactions
  Interaction Bpp  { type = zeeman;   spins = e,h; field = "0 0 -0.2005"; }
  Interaction HF_e { type = hyperfine; group1 = e; group2 = n1; prefactor = 1e-3; } // uses n1 tensor
  Interaction HF_h { type = hyperfine; group1 = h; group2 = n2; prefactor = 1e-3; } // uses n2 tensor
  //Interaction Jpp  { type = doublespin; group1 = e; group2 = h; tensor = isotropic("0.0"); prefactor = 1e-3; }

  // ---------- PP electronic ⊗ nuclear states (uu/ud/du/dd) ----------
  // singlet and triplet sublevels (explicit); nuclear tags preserved
  // uu
  State PP_S_uu  { spins(e,h)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1)=| 1/2>; spin(n2)=| 1/2>; }
  State PP_T0_uu { spins(e,h)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1)=| 1/2>; spin(n2)=| 1/2>; }
  State PP_Tp_uu { spin(e)=| 1/2>; spin(h)=| 1/2>;        spin(n1)=| 1/2>; spin(n2)=| 1/2>; }
  State PP_Tm_uu { spin(e)=|-1/2>; spin(h)=|-1/2>;        spin(n1)=| 1/2>; spin(n2)=| 1/2>; }
  // ud
  State PP_S_ud  { spins(e,h)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1)=| 1/2>; spin(n2)=|-1/2>; }
  State PP_T0_ud { spins(e,h)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1)=| 1/2>; spin(n2)=|-1/2>; }
  State PP_Tp_ud { spin(e)=| 1/2>; spin(h)=| 1/2>;        spin(n1)=| 1/2>; spin(n2)=|-1/2>; }
  State PP_Tm_ud { spin(e)=|-1/2>; spin(h)=|-1/2>;        spin(n1)=| 1/2>; spin(n2)=|-1/2>; }
  // du
  State PP_S_du  { spins(e,h)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1)=|-1/2>; spin(n2)=| 1/2>; }
  State PP_T0_du { spins(e,h)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1)=|-1/2>; spin(n2)=| 1/2>; }
  State PP_Tp_du { spin(e)=| 1/2>; spin(h)=| 1/2>;        spin(n1)=|-1/2>; spin(n2)=| 1/2>; }
  State PP_Tm_du { spin(e)=|-1/2>; spin(h)=|-1/2>;        spin(n1)=|-1/2>; spin(n2)=| 1/2>; }
  // dd
  State PP_S_dd  { spins(e,h)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1)=|-1/2>; spin(n2)=|-1/2>; }
  State PP_T0_dd { spins(e,h)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1)=|-1/2>; spin(n2)=|-1/2>; }
  State PP_Tp_dd { spin(e)=| 1/2>; spin(h)=| 1/2>;        spin(n1)=|-1/2>; spin(n2)=|-1/2>; }
  State PP_Tm_dd { spin(e)=|-1/2>; spin(h)=|-1/2>;        spin(n1)=|-1/2>; spin(n2)=|-1/2>; }

  // Initial PP population: F-pair (equal mixture)
  Properties prop { InitialState =  PP_S_uu, PP_T0_uu, PP_Tp_uu, PP_Tm_uu, PP_S_ud, PP_T0_ud, PP_Tp_ud, PP_Tm_ud, PP_S_du, PP_T0_du, PP_Tp_du, PP_Tm_du, PP_S_dd, PP_T0_dd, PP_Tp_dd, PP_Tm_dd;}

  // ---------- PP → X formation (spin-selective; preserve nuclear tag) ----------
  Transition formS_uu  { source = PP_S_uu;  targetsystem = X; targetstate = X_S1_uu; rate = 0.060; }
  Transition formS_ud  { source = PP_S_ud;  targetsystem = X; targetstate = X_S1_ud; rate = 0.060; }
  Transition formS_du  { source = PP_S_du;  targetsystem = X; targetstate = X_S1_du; rate = 0.060; }
  Transition formS_dd  { source = PP_S_dd;  targetsystem = X; targetstate = X_S1_dd; rate = 0.060; }

  Transition formT0_uu { source = PP_T0_uu; targetsystem = X; targetstate = X_T0_uu; rate = 0.020; }
  Transition formT0_ud { source = PP_T0_ud; targetsystem = X; targetstate = X_T0_ud; rate = 0.020; }
  Transition formT0_du { source = PP_T0_du; targetsystem = X; targetstate = X_T0_du; rate = 0.020; }
  Transition formT0_dd { source = PP_T0_dd; targetsystem = X; targetstate = X_T0_dd; rate = 0.020; }

  Transition formTp_uu { source = PP_Tp_uu; targetsystem = X; targetstate = X_Tp_uu; rate = 0.020; }
  Transition formTp_ud { source = PP_Tp_ud; targetsystem = X; targetstate = X_Tp_ud; rate = 0.020; }
  Transition formTp_du { source = PP_Tp_du; targetsystem = X; targetstate = X_Tp_du; rate = 0.020; }
  Transition formTp_dd { source = PP_Tp_dd; targetsystem = X; targetstate = X_Tp_dd; rate = 0.020; }

  Transition formTm_uu { source = PP_Tm_uu; targetsystem = X; targetstate = X_Tm_uu; rate = 0.020; }
  Transition formTm_ud { source = PP_Tm_ud; targetsystem = X; targetstate = X_Tm_ud; rate = 0.020; }
  Transition formTm_du { source = PP_Tm_du; targetsystem = X; targetstate = X_Tm_du; rate = 0.020; }
  Transition formTm_dd { source = PP_Tm_dd; targetsystem = X; targetstate = X_Tm_dd; rate = 0.020; }

  // Operator lindbladrelaxation {type=relaxationLindblad; rate1=0.001; rate2=0.001; rate3=0.001; spins=e,h;}

}

// ---------------------- X: exciton manifold (S1, T0, Tp, Tm) ----------------------
SpinSystem X
{
  // bound pair (no Δg in X)
  Spin eX  { type = electron; spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 2.0023"); }
  Spin hX  { type = electron; spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 2.0023"); }
  Spin n1X { spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 1");} 
  Spin n2X { spin = 1/2; tensor = matrix("0 0 0; 0 0 0; 0 0 1");}

  // Zeeman (swept in lockstep with PP)
  Interaction Bx  { type = zeeman; spins = eX,hX; field = "0 0 -0.2005"; }

  // Triplet ZFS proxy — principal values in mT (choose per emitter)
  // 3CzClIPN:  -28  -52   80   (HF knee ~80 mT)
  // 4CzIPN:    -15  -30   45   (HF knee ~60 mT)   ← use these three numbers for 4CzIPN
  // 5CzBN:     -7.5 -15  22.5  (HF knee ~30 mT)   ← use these three numbers for 5CzBN
  //Interaction Dex { type = doublespin; group1 = eX; group2 = hX; tensor = anisotropic("-28 -52 80"); prefactor = 2.0023e-3; ignoretensors=true;}
  //  Interaction Dex { type = doublespin; group1 = eX; group2 = hX; tensor = matrix("0.000000 0.000000 25.980762; 0.000000 -30.000000 0.000000; 25.980762 0.000000 30.000000"); prefactor = 2.0023e-3; ignoretensors = true;}
  //  Interaction Dex { type = doublespin; group1 = eX; group2 = hX; tensor = matrix("-1.000000 0.000000 46.765372; 0.000000 -52.000000 0.000000; 46.765372 0.000000 53.000000"); prefactor = 2.0023e-3; ignoretensors = true;}
  Interaction Dex { type = doublespin; group1 = eX; group2 = hX; tensor = matrix("0.000000 0.000000 51.961524; 0.000000 -60.000000 0.000000; 51.961524 0.000000 60.000000"); prefactor = 2.0023e-3; ignoretensors = true;}

  // no residual exchange inside X
  // Interaction Jx  { type = doublespin; group1 = eX; group2 = hX; tensor = isotropic("0.0"); prefactor = 1e-3; }

  // ---------- X states (declare all; names must match PP→X targets) ----------
  // uu
  State X_S1_uu { spins(eX,hX)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1X)=| 1/2>; spin(n2X)=| 1/2>; }
  State X_T0_uu { spins(eX,hX)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1X)=| 1/2>; spin(n2X)=| 1/2>; }
  State X_Tp_uu { spin(eX)=| 1/2>; spin(hX)=| 1/2>;        spin(n1X)=| 1/2>; spin(n2X)=| 1/2>; }
  State X_Tm_uu { spin(eX)=|-1/2>; spin(hX)=|-1/2>;        spin(n1X)=| 1/2>; spin(n2X)=| 1/2>; }
  // ud
  State X_S1_ud { spins(eX,hX)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1X)=| 1/2>; spin(n2X)=|-1/2>; }
  State X_T0_ud { spins(eX,hX)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1X)=| 1/2>; spin(n2X)=|-1/2>; }
  State X_Tp_ud { spin(eX)=| 1/2>; spin(hX)=| 1/2>;        spin(n1X)=| 1/2>; spin(n2X)=|-1/2>; }
  State X_Tm_ud { spin(eX)=|-1/2>; spin(hX)=|-1/2>;        spin(n1X)=| 1/2>; spin(n2X)=|-1/2>; }
  // du
  State X_S1_du { spins(eX,hX)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1X)=|-1/2>; spin(n2X)=| 1/2>; }
  State X_T0_du { spins(eX,hX)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1X)=|-1/2>; spin(n2X)=| 1/2>; }
  State X_Tp_du { spin(eX)=| 1/2>; spin(hX)=| 1/2>;        spin(n1X)=|-1/2>; spin(n2X)=| 1/2>; }
  State X_Tm_du { spin(eX)=|-1/2>; spin(hX)=|-1/2>;        spin(n1X)=|-1/2>; spin(n2X)=| 1/2>; }
  // dd
  State X_S1_dd { spins(eX,hX)=| 1/2,-1/2> - |-1/2, 1/2>;  spin(n1X)=|-1/2>; spin(n2X)=|-1/2>; }
  State X_T0_dd { spins(eX,hX)=| 1/2,-1/2> + |-1/2, 1/2>;  spin(n1X)=|-1/2>; spin(n2X)=|-1/2>; }
  State X_Tp_dd { spin(eX)=| 1/2>; spin(hX)=| 1/2>;        spin(n1X)=|-1/2>; spin(n2X)=|-1/2>; }
  State X_Tm_dd { spin(eX)=|-1/2>; spin(hX)=|-1/2>;        spin(n1X)=|-1/2>; spin(n2X)=|-1/2>; }

  // ---------- RISC: T → S1 (field-independent) ----------
  //Transition RISC0_uu { source = X_T0_uu; targetstate = X_S1_uu; rate = 0.04; }
  //Transition RISCp_uu { source = X_Tp_uu; targetstate = X_S1_uu; rate = 0.04; }
  //Transition RISCm_uu { source = X_Tm_uu; targetstate = X_S1_uu; rate = 0.04; }
  //Transition RISC0_ud { source = X_T0_ud; targetstate = X_S1_ud; rate = 0.04; }
  //Transition RISCp_ud { source = X_Tp_ud; targetstate = X_S1_ud; rate = 0.04; }
  //Transition RISCm_ud { source = X_Tm_ud; targetstate = X_S1_ud; rate = 0.04; }
  //Transition RISC0_du { source = X_T0_du; targetstate = X_S1_du; rate = 0.04; }
  //Transition RISCp_du { source = X_Tp_du; targetstate = X_S1_du; rate = 0.04; }
  //Transition RISCm_du { source = X_Tm_du; targetstate = X_S1_du; rate = 0.04; }
  //Transition RISC0_dd { source = X_T0_dd; targetstate = X_S1_dd; rate = 0.04; }
  //Transition RISCp_dd { source = X_Tp_dd; targetstate = X_S1_dd; rate = 0.04; }
  //Transition RISCm_dd { source = X_Tm_dd; targetstate = X_S1_dd; rate = 0.04; }

  // ---------- TCA: Triplet-charge annihilation (non-radiative loss) ----------
  // 3CzClIPN / 4CzIPN: TCA0 ≫ TCA±  → positive HF branch
  // 5CzBN: keep TCA very small (we want TTA-dominated shape)
  Transition TCA0_uu { source = X_T0_uu; targetsystem = SINK; targetstate = NR; rate = 0.020; }
  Transition TCAp_uu { source = X_Tp_uu; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCAm_uu { source = X_Tm_uu; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCA0_ud { source = X_T0_ud; targetsystem = SINK; targetstate = NR; rate = 0.020; }
  Transition TCAp_ud { source = X_Tp_ud; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCAm_ud { source = X_Tm_ud; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCA0_du { source = X_T0_du; targetsystem = SINK; targetstate = NR; rate = 0.020; }
  Transition TCAp_du { source = X_Tp_du; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCAm_du { source = X_Tm_du; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCA0_dd { source = X_T0_dd; targetsystem = SINK; targetstate = NR; rate = 0.020; }
  Transition TCAp_dd { source = X_Tp_dd; targetsystem = SINK; targetstate = NR; rate = 0.0005; }
  Transition TCAm_dd { source = X_Tm_dd; targetsystem = SINK; targetstate = NR; rate = 0.0005; }

  // ---------- TTA: Triplet-triplet annihilation → radiative yield ----------
  // (keep tiny for 3CzClIPN/4CzIPN; make it dominant for 5CzBN)
  // 3CzClIPN defaults below; for 5CzBN set TTA0=0.080 and TTA±=0.005; for 4CzIPN set TTA0=0.008 and TTA±=0.002
  Transition TTA0_uu { source = X_T0_uu; targetsystem = SINK; targetstate = RAD; rate = 0.005; }
  Transition TTAp_uu { source = X_Tp_uu; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTAm_uu { source = X_Tm_uu; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTA0_ud { source = X_T0_ud; targetsystem = SINK; targetstate = RAD; rate = 0.005; }
  Transition TTAp_ud { source = X_Tp_ud; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTAm_ud { source = X_Tm_ud; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTA0_du { source = X_T0_du; targetsystem = SINK; targetstate = RAD; rate = 0.005; }
  Transition TTAp_du { source = X_Tp_du; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTAm_du { source = X_Tm_du; targetsystem = SINK; targetstate = RAD; rate = 0.001; } // (typo guard: ensure names match)
  Transition TTA0_dd { source = X_T0_dd; targetsystem = SINK; targetstate = RAD; rate = 0.005; }
  Transition TTAp_dd { source = X_Tp_dd; targetsystem = SINK; targetstate = RAD; rate = 0.001; }
  Transition TTAm_dd { source = X_Tm_dd; targetsystem = SINK; targetstate = RAD; rate = 0.001; }

  // ---------- S1 decay ----------
  Transition kRad_uu  { source = X_S1_uu; targetsystem = SINK; targetstate = RAD; rate = 0.010; }
  Transition kNR_S_uu { source = X_S1_uu; targetsystem = SINK; targetstate = NR;  rate = 0.002; }
  Transition kRad_ud  { source = X_S1_ud; targetsystem = SINK; targetstate = RAD; rate = 0.010; }
  Transition kNR_S_ud { source = X_S1_ud; targetsystem = SINK; targetstate = NR;  rate = 0.002; }
  Transition kRad_du  { source = X_S1_du; targetsystem = SINK; targetstate = RAD; rate = 0.010; }
  Transition kNR_S_du { source = X_S1_du; targetsystem = SINK; targetstate = NR;  rate = 0.002; }
  Transition kRad_dd  { source = X_S1_dd; targetsystem = SINK; targetstate = RAD; rate = 0.010; }
  Transition kNR_S_dd { source = X_S1_dd; targetsystem = SINK; targetstate = NR;  rate = 0.002; }

  //Operator lindbladrelaxation {type=relaxationLindblad; rate1=0.001; rate2=0.001; rate3=0.001; spins=eX,hX;}
}

// ---------------------- SINK: radiative / non-radiative counters ----------------------
SpinSystem SINK
{
  // two inert spins to host distinct sink states (naming convenience)
  Spin s1 { spin = 1/2; } 
  Spin s2 { spin = 1/2; }
  Spin s3 { spin = 1/2; }
  Spin s4 { spin = 1/2; }
  State RAD { spin(s1)=| 1/2>; spin(s2)=| 1/2>; }
  State NR  { spin(s1)=|-1/2>; spin(s2)=|-1/2>; }
}

// ---------------------- Settings (sweep & outputs) ----------------------
Settings
{
  Settings general { steps = 1; }

  // sweep both PP and X Zeeman fields by +0.001 T each step (x,y,z in Tesla)
  Action sweep_pp { type = addvector; vector = PP.Bpp.field; direction = "0 0 1"; value = 0.001; }
  Action sweep_x  { type = addvector; vector = X.Bx.field;   direction = "0 0 1"; value = 0.001; }

  // export field components so your script can read B (pp.bpp.field.z)
  Output Bfield { type = xyz; vector = PP.Bpp.field; }
}

// ---------------------- Run (finite-time propagation; sinks accumulate) ----------------------
Run
{
  Task TE
  {
    type       = MultiStaticSS-TimeEvolution;
    logfile    = "run_3czclpn.log";
    datafile   = "run_3czclpn.dat";

    // stable per-step Δyields; pick either (A) or (B):
    // (A) "fast": 
    //TimeStep   = 1;     // ns
    //TotalTime  = 3000;     // ns

    // (B) "fine":
    timestep   = 0.1;   // ns
    totaltime  = 8000;    // ns
    //maximumtimestep = 0.005;

    settings   = general;
  }
}


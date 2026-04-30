#!/usr/bin/env wolframscript
(* Times the Mathematica BenderWu.m implementation across the same matrix as
   benchmark/cases.jl and writes benchmark/results_mma.json.

   Usage:
     wolframscript -file benchmark/run_mathematica.wls [--quick]

   On macOS without wolframscript on PATH, invoke via the full path:
     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file ... *)

Get["~/Library/Wolfram/Applications/BenderWu.m"];

quick = MemberQ[$ScriptCommandLine, "--quick"];

(* Mirror benchmark/cases.jl. Keep coefficients exact (rationals/integers) so
   the "exact" mode runs symbolically. The Float64 mode wraps the potential in
   N[..., MachinePrecision] before calling BenderWu. *)
potentials = {
  <|"name" -> "quartic",      "V" -> x^2/2 + x^4|>,
  <|"name" -> "sextic",       "V" -> x^2/2 + x^6|>,
  <|"name" -> "octic",        "V" -> x^2/2 + x^8|>,
  <|"name" -> "mixed_parity", "V" -> x^2/2 + x^3 + x^4|>
};

nus            = If[quick, {0},     {0, 1, 5}];
ordersFloat    = If[quick, {10},    {10, 20, 30, 40, 50}];
ordersRational = If[quick, {5},     {5, 10, 15, 20, 25}];

(* Run BenderWu and return the energy ε list as {e0, ε_2, ε_4, ...}.
   Output -> "Energy", OutputStyle -> "Array" returns Drop[epsilon, lmin/2+1]
   per BenderWu.m line 247. lmin defaults to 0 so we get a flat list. *)
runOnce[V_, nu_, n_, prec_] := Module[{Vp},
  Vp = If[prec === "exact", V, N[V, MachinePrecision]];
  (* Monitor -> False suppresses BenderWu's per-order Dynamic progress
     panel, which is meaningless under wolframscript and floods stderr.
     Warnings -> False silences the "energy expansion contains zeros"
     notes that fire for sextic/octic where odd orders vanish. *)
  BenderWu[Vp, x, nu, n,
    Output -> "Energy", OutputStyle -> "Array",
    Monitor -> False, Warnings -> False]
];

(* RepeatedTiming chooses repetitions automatically and returns the average.
   For very fast cells it averages many runs; for slow ones it runs once.
   Returns {avgSeconds, result}. *)
timeOnce[V_, nu_, n_, prec_] := RepeatedTiming[runOnce[V, nu, n, prec]];

(* Serialise an exact value as Mathematica InputForm so it round-trips. For
   machine-precision floats we want a high-precision decimal so the Julia
   side can compare without re-parsing Mathematica syntax. *)
serialiseEpsilon[v_, "exact"] := ToString[v, InputForm];
serialiseEpsilon[v_, _] := ToString[CForm[N[v, 20]]];

results = {};
i = 0;
total = Length[potentials] * Length[nus] *
        (Length[ordersFloat] + Length[ordersRational]);

Do[
  Module[{Vexpr, name, nu, N0, prec, orders, t, eps, dt, sample},
    name  = pot["name"];
    Vexpr = pot["V"];
    Do[
      Do[
        {prec, orders} = modeSpec;
        Do[
          i++;
          {dt, eps} = timeOnce[Vexpr, nu, N0, prec];
          (* RepeatedTiming returns wall-clock seconds averaged over repeats;
             we record it as both "median" (proxy) and "min" (same value). *)
          AppendTo[results, <|
            "potential" -> name,
            "nu" -> nu,
            "N" -> N0,
            "mode" -> If[prec === "exact", "rational", "float64"],
            "time_ns_median" -> dt * 10^9,
            "time_ns_min"    -> dt * 10^9,
            "epsilons" -> Map[serialiseEpsilon[#, prec] &, eps]
          |>];
          Print["[", i, "/", total, "] ", name, " nu=", nu, " N=", N0,
                " mode=", If[prec === "exact", "rational", "float64"],
                " ... ", ToString[NumberForm[dt * 1000.0, {7, 3}]], " ms"];
          , {N0, orders}]
        , {modeSpec, {{"machine", ordersFloat}, {"exact", ordersRational}}}]
      , {nu, nus}]
  ]
  , {pot, potentials}];

machine = <|
  "mathematica_version" -> $Version,
  "kernel_id" -> $MachineName,
  "cpu" -> $ProcessorType,
  "os" -> $OperatingSystem
|>;

out = <|"machine" -> machine, "results" -> results|>;

(* Write JSON. Use Compact->2 for some pretty-printing without huge files. *)
outPath = FileNameJoin[{DirectoryName[$InputFileName], "results_mma.json"}];
Export[outPath, out, "JSON"];
Print["Wrote ", outPath];

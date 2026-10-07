#!/usr/bin/env wolframscript
(* Times the Mathematica BenderWu.m implementation in exact-rational mode
   across the same matrix as benchmark/cases.jl and writes
   benchmark/results_mma.json.

   Usage:
     wolframscript -file benchmark/run_mathematica.wls [--quick]

   On macOS without wolframscript on PATH, invoke via the full path:
     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file ... *)

Get["~/Library/Wolfram/Applications/BenderWu.m"];

quick = MemberQ[$ScriptCommandLine, "--quick"];

(* Mirror benchmark/cases.jl. Coefficients stay exact (rationals/integers) so
   the computation runs symbolically. Orders are N, the highest power of g;
   BenderWu counts powers of g^2, so it is called with N/2. *)
potentials = {
  <|"name" -> "quartic",      "V" -> x^2/2 + x^4|>,
  <|"name" -> "mixed_parity", "V" -> x^2/2 + x^3 + x^4|>,
  <|"name" -> "sextic",       "V" -> x^2/2 + x^6|>,
  <|"name" -> "octic",        "V" -> x^2/2 + x^8|>
};

nus    = If[quick, {0}, {0, 1, 5}];
orders = If[quick, {10}, {10, 20, 30, 40, 50}];

(* Run BenderWu up to g^n0 and return the n0/2 + 1 energy coefficients
   {ε_0, ε_2, ..., ε_n0}. BenderWu's order argument counts powers of g^2.
   Output -> "Energy", OutputStyle -> "Array" returns Drop[epsilon, lmin/2+1]
   per BenderWu.m line 247. lmin defaults to 0 so we get a flat list.
   Monitor -> False suppresses BenderWu's per-order Dynamic progress panel.
   Warnings -> False silences the "expansion contains zeros" notes for
   sextic/octic where odd orders vanish. *)
runOnce[V_, nu_, n0_] := BenderWu[V, x, nu, n0/2,
  Output -> "Energy", OutputStyle -> "Array",
  Monitor -> False, Warnings -> False];

(* Mirror the BenchmarkTools settings in run_julia.jl: after one warm-up
   call, collect AbsoluteTiming samples until the 5 s budget is spent or 50
   samples are taken, with at least one sample. Returns
   {sampleSeconds, result}. *)
timeCase[V_, nu_, n0_] := Module[{eps, times = {}, start},
  eps = runOnce[V, nu, n0];
  start = AbsoluteTime[];
  While[times === {} || (Length[times] < 50 && AbsoluteTime[] - start < 5),
    AppendTo[times, First[AbsoluteTiming[runOnce[V, nu, n0]]]]];
  {times, eps}];

results = {};
i = 0;
total = Length[potentials] * Length[nus] * Length[orders];

Do[
  Module[{Vexpr, name, nu, N0, times, eps, dt},
    name  = pot["name"];
    Vexpr = pot["V"];
    Do[
      Do[
        i++;
        {times, eps} = timeCase[Vexpr, nu, N0];
        dt = Median[times];
        AppendTo[results, <|
          "potential" -> name,
          "nu" -> nu,
          "N" -> N0,
          "time_ns_median" -> dt * 10^9,
          "time_ns_min"    -> Min[times] * 10^9,
          "samples"        -> Length[times],
          "epsilons" -> Map[ToString[#, InputForm] &, eps]
        |>];
        Print["[", i, "/", total, "] ", name, " nu=", nu, " N=", N0,
              " ... ", ToString[NumberForm[dt * 1000.0, {7, 3}]], " ms (median)"];
        , {N0, orders}]
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

outPath = FileNameJoin[{DirectoryName[$InputFileName], "results_mma.json"}];
Export[outPath, out, "JSON"];
Print["Wrote ", outPath];

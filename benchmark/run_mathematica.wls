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
   the computation runs symbolically. *)
potentials = {
  <|"name" -> "quartic",      "V" -> x^2/2 + x^4|>,
  <|"name" -> "sextic",       "V" -> x^2/2 + x^6|>,
  <|"name" -> "octic",        "V" -> x^2/2 + x^8|>,
  <|"name" -> "mixed_parity", "V" -> x^2/2 + x^3 + x^4|>
};

nus    = If[quick, {0}, {0, 1, 5}];
orders = If[quick, {5}, {5, 10, 15, 20, 25}];

(* Run BenderWu and return the energy ε list as {e0, ε_2, ε_4, ...}.
   Output -> "Energy", OutputStyle -> "Array" returns Drop[epsilon, lmin/2+1]
   per BenderWu.m line 247. lmin defaults to 0 so we get a flat list.
   Monitor -> False suppresses BenderWu's per-order Dynamic progress panel.
   Warnings -> False silences the "expansion contains zeros" notes for
   sextic/octic where odd orders vanish. *)
runOnce[V_, nu_, n_] := BenderWu[V, x, nu, n,
  Output -> "Energy", OutputStyle -> "Array",
  Monitor -> False, Warnings -> False];

(* RepeatedTiming chooses repetitions automatically and returns the average.
   For very fast cells it averages many runs; for slow ones it runs once.
   Returns {avgSeconds, result}. *)
timeOnce[V_, nu_, n_] := RepeatedTiming[runOnce[V, nu, n]];

results = {};
i = 0;
total = Length[potentials] * Length[nus] * Length[orders];

Do[
  Module[{Vexpr, name, nu, N0, t, eps, dt},
    name  = pot["name"];
    Vexpr = pot["V"];
    Do[
      Do[
        i++;
        {dt, eps} = timeOnce[Vexpr, nu, N0];
        AppendTo[results, <|
          "potential" -> name,
          "nu" -> nu,
          "N" -> N0,
          "time_ns_median" -> dt * 10^9,
          "time_ns_min"    -> dt * 10^9,
          "epsilons" -> Map[ToString[#, InputForm] &, eps]
        |>];
        Print["[", i, "/", total, "] ", name, " nu=", nu, " N=", N0,
              " ... ", ToString[NumberForm[dt * 1000.0, {7, 3}]], " ms"];
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

---
layout: default
title: Options
---

The Options section is enclosed between

    <BeginOptions>
    ...
    <EndOptions>

The following option keys and values are known.

|Option Key|Option Value|Action|
|----------|------------|------|
|Accuracy|VeryTight|Converge all parts of the SCF to a very tight accuracy|
|BasisSets|6-31G\*|Use a 6-31G\*\* basis set|
|Charge|Integer|The total charge of the system|
|DebugAll|MaxDebug|Set the debugging to maximum|
|Geometry|InAngstrom|The coordinates are in units of Angstrom|
|Grad|Optimize|Optimize the configuration|
|Guess|Superpos|Make a superposition guess|
|ModelChem|B3LYP|Use the B3LYP exchange-correlation functional|
|Multiplicity|Integer|The multiplicity of the system|
|OutPut|XYZ|The format of the output configuration file|
|SCFConvergence|DIIS|Use the DIIS algorithm to speedup SCF convergence|
|SCFMethod|RH|Use Roothan-Hall|

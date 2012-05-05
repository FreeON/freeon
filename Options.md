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
|Grad|Optimize|Optimize the configuration|
|Charge|Integer|The total charge of the system|
|Multiplicity|Integer|The multiplicity of the system|
|Guess|Superpos|Make a superposition guess|
|OutPut|XYZ|The format of the output configuration file|
|DebugAll|MaxDebug|Set the debugging to maximum|
|SCFMethod|RH|Use Roothan-Hall|
|BasisSets|6-31G\*|Use a 6-31G\*\* basis set|
|ModelChem|B3LYP|Use the B3LYP exchange-correlation functional|
|Accuracy|VeryTight|Converge all parts of the SCF to a very tight accuracy|
|SCFConvergence|DIIS|Use the DIIS algorithm to speedup SCF convergence|
|Geometry|InAngstrom|The coordinates are in units of Angstrom|



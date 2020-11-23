# reference from https://www.bioinfor.com/amino-acid/
const RESNAME2MASS = Dict{String, Float32}(
"ALA" =>  71.03711,
"ARG" => 156.10111,
"ASN" => 114.04293,
"ASP" => 115.02694,
"CYS" => 103.00919,
"GLU" => 129.04259,
"GLN" => 128.05858,
"GLY" =>  57.02146,
"HIS" => 137.05895,
"ILE" => 113.08406,
"LEU" => 113.08406,
"LYS" => 128.09496,
"MET" => 131.04049,
"PHE" => 147.06841,
"PRO" =>  97.05276,
"SER" =>  87.03203,
"THR" => 101.04768,
"TRP" => 186.07931,
"TYR" => 163.06333,
"VAL" =>  99.06841,
"HOH" =>  18.01528
                                    )

const ATOMNAME2MASS = Dict{String, Float32}(
"H"  => 1.00798,
"C"  => 12.0106,
"N"  => 14.0069,
"O"  => 15.9994,
"P"  => 30.9738,
"S"  => 32.0680,
"CL" => 35.4520,
"NA" => 22.9898
                                     )

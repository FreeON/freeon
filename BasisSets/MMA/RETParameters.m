(* H parameters *)
 
b[s,H]=-0.425;bp[s,H]=0.928;a[s,H]=0.3243;ap[s,H]=-3.692;
 
(* O parameters *)
 
b[s,O]=-0.5063;bp[s,O]=1.2788;a[s,O]=0.543;ap[s,O]=-2.7962;
b[p,O]=-0.4413;bp[p,O]=0.9128;a[p,O]=0.6088;ap[p,O]=-3.2801;

(* C parameters *)

(* old parameters ....
 
b[s,C]=-0.5120;bp[s,C]=1.2953;a[s,C]=0.4782;ap[s,C]=-3.3667;
b[p,C]=-0.4361;bp[p,C]=0.8857;a[p,C]=0.5278;ap[p,C]=-3.7305;

*)


dilate=8/10;
b[s,C]=-0.5120;bp[s,C]=1.2953;a[s,C]=0.4782;ap[s,C]=-3.3667*dilate;
b[p,C]=-0.4361;bp[p,C]=0.8857;a[p,C]=0.5278;ap[p,C]=-3.7305*dilate;

.. _bands2xmgr:

bands2xmgr
==========

Utility for visualizing SIESTA bands and density-of-states results.

The function bands2xmgr reads "systemlabel.bands" and
"systemlabel.DOS" of SIESTA and writes "systemlabel.agr" for XMGR/GRACE
and (optionally) an eps file "systemlabel.eps".

Usage :
  bands2xmgr systemlabel [flags]

where the optional keyword flags are:

    -S   : Sort E(k) for smoother band lines,
    -O   : Output (sorted) band structure to specified file ,
    -R   : Sort in reverse order,
    -fs  : Set the flip sensitivity factor (0.<fs<=1.) where fs=1 corresponds to  maximum sensitivity (many band flips),
    -f   : Set first band index in the sorting,
    -l   : Set last band index in the sorting.
    -P   : Print eps file

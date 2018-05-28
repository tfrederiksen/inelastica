Minimal NEB example
===================

NH3 molecule inversion in 10 Ang^3 box.
Start calculation with:

      NEB -n 9 -p 8 -c'2, 3,-0.476,0.824,0, 1,-0.476,0.824,0'  L R

Since we have a small molecule we constrain the geometry (to fix translational and rotational symmetries).
First we keep one of the hydrogens fixed `'2,'` (where 2 is the number of the hydrogen atom).
Then a second hydrogen is fixed along the vector `'3,-0.476,0.824,0,'` (3 is the number of the hydrogen
and after comes the vector, see the `STRUCT.fdf` to see why I choose this vector)
and finally, the nitrogen is allowed to move in the plane perpendicular to the vector, `'1,-0.476,0.824,0'` (1 is the nitrogen number).


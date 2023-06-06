/*
Rectangular mesh
*/

lc = 2; // mesh precision
Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 100, 0, lc};
Point(4) = {0, 100, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Physical Line("N") = {1, 3};
Physical Line("D") = {2, 4};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface("S") = {1};

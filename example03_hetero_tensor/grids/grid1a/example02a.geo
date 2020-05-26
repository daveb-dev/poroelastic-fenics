Mesh.MshFileVersion = 2.8;
//+
l = 10;
b = 10;
h = 2;
//+
radius = 0.6096;
size = 0.6;
refine2=2;
//+
Point(1) = {0, 0, 0, refine2};
//+
Point(2) = {radius, 0, 0, size};
//+
Point(3) = {l, 0, 0, refine2};
//+
Point(4) = {l, b, 0, refine2};
//+
Point(5) = {0, b, 0, refine2};
//+
Point(6) = {0, radius, 0, size};


//+
Circle(1) = {6, 1, 2};
//+
Line(2) = {6, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {4, 3};
//+
Line(5) = {3, 2};
//+
Curve Loop(1) = {2, 3, 4, 5, -1};
//+
Plane Surface(1) = {1};


//+
Extrude {0, 0, h} {
  Surface{1}; Recombine;
}

//+
Physical Point(1) = {7};
//+
Physical Point(2) = {20};
//+
Physical Point(3) = {16};
//+
Physical Point(4) = {12};
//+
Physical Point(5) = {8};
//+
Physical Point(6) = {5};
//+
Physical Point(7) = {6};
//+
Physical Point(8) = {2};
//+
Physical Point(9) = {3};
//+
Physical Point(10) = {4};
//+
Physical Curve(11) = {5};
//+
Physical Curve(12) = {10};
//+
Physical Curve(13) = {22};
//+
Physical Curve(14) = {9};
//+
Physical Curve(15) = {4};
//+
Physical Curve(16) = {18};
//+
Physical Curve(17) = {3};
//+
Physical Curve(18) = {8};
//+
Physical Curve(19) = {7};
//+
Physical Curve(20) = {14};
//+
Physical Curve(21) = {2};
//+
Physical Curve(22) = {13};
//+
Physical Curve(23) = {26};
//+
Physical Curve(24) = {1};
//+
Physical Curve(25) = {11};
//+
Physical Surface(26) = {32};
//+
Physical Surface(27) = {19};
//+
Physical Surface(28) = {15};
//+
Physical Surface(29) = {1};
//+
Physical Surface(30) = {27};
//+
Physical Surface(31) = {23};
//+
Physical Volume(32) = {1};


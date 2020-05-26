//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 0.2, 0, 1.0};
//+
Point(4) = {0, 0.5, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {1, 0.2, 0, 1.0};
//+
Point(7) = {1, 0.5, 0, 1.0};
//+
Point(8) = {1, 1, 0, 1.0};
//+
Point(9) = {0, 0.325, 0, 1.0};
//+
Point(10) = {0, 0.375, 0, 1.0};

//+
Line(1) = {3, 6};
//+
Line(2) = {6, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 3};
//+
Line(5) = {4, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {3, 9};
//+
Line(8) = {9, 10};
//+
Line(9) = {10, 4};
//+
Line(10) = {5, 8};
//+
Line(11) = {8, 7};
//+
Line(12) = {4, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, -6, -5, -9, -8, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, -11, -10, -12};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {3, 3, 10} = 12 Using Progression 1.1;
//+
Transfinite Curve {1, 5} = 20 Using Progression 1.1;
//+
Transfinite Curve {9, 7} = 5 Using Progression 1;
//+
Transfinite Curve {8} = 5 Using Progression 1;
//+
Transfinite Curve {12} = 10 Using Progression 1.1;
//+
Transfinite Curve {11} = 7 Using Progression 0.9;
//+
Transfinite Curve {4} = 5 Using Progression 0.9;
//+
Transfinite Curve {2} = 5 Using Progression 0.9;
//+
Transfinite Curve {6} = 5 Using Progression 1;
//+
Physical Curve("well") = {8};
//+
Physical Curve("bottom") = {3};
//+
Physical Curve("right") = {2, 6, 11};
//+
Physical Curve("top") = {10};
//+
Physical Curve("left") = {12, 9, 7, 4};
//+
Physical Surface("top") = {3};
//+
Physical Surface("reservoir") = {2};
//+
Physical Surface("bottom") = {1};

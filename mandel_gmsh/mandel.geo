//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Line("left") = {4};
//+
Physical Line("bottom") = {1};
//+
Physical Line("right") = {2};
//+
Physical Line("top") = {3};
//+
Transfinite Line {3} = 20 Using Progression 1.1;
//+
Transfinite Line {1} = 20 Using Progression 0.9;
//+
Transfinite Line {2} = 40 Using Progression 1;
//+
Transfinite Line {4} = 20 Using Progression 1;
//+
Plane Surface(2) = {1};
//+
Physical Surface("surface1") = {1};




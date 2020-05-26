//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.05, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, 1, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {0, 0.05, 0, 1.0};


//+
Circle(1) = {6, 1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Curve Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.2} {
  Curve{4}; Curve{3}; Curve{2}; Surface{1}; Curve{5}; Curve{1}; 
}
//+
Extrude {0, 0, 0.4} {
  Surface{44}; // Layers{3}; Recombine;
}
//+
Extrude {0, 0, 0.4} {
  Surface{71}; //Layers{3}; Recombine;
}

//+
Physical Surface("top") = {98};
//+
Physical Surface("bottom") = {1};
//+
Physical Surface("front") = {17,66,93};
//+
Physical Surface("right") = {13,70,97};
//+
Physical Surface("back") = {9,54,81};
//+
Physical Surface("left") = {31,58,85};
//+
Physical Surface("well") = {62};
//+
Physical Surface("well_noflow") = {35,89};
//+
Physical Volume("top_layer") = {1};
//+
Physical Volume("middle_layer") = {2};
//+
Physical Volume("bottom_layer") = {3};


//+
Transfinite Curve {2, 14, 49, 76} = 10 Using Progression 1.1;
//+
Transfinite Curve {5, 20, 47, 74} = 10 Using Progression 0.9;
//+
Transfinite Curve {3, 10, 50, 77} = 10 Using Progression 1.0;
//+
Transfinite Curve {4, 6, 46, 73} = 10 Using Progression 0.9;
//+
Transfinite Curve {7, 11, 15, 30, 8} = 5 Using Progression 0.9;
//+
Transfinite Curve {53, 57, 61, 65, 52} = 10 Using Progression 1;
//+
Transfinite Curve {80, 84, 88, 92, 79} = 5 Using Progression 1.1;
//+
Transfinite Curve {1, 21, 48, 75} = 15 Using Progression 1;
//+
Physical Surface(12) = {89, 62, 35};

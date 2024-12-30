R = 0.3;//设定圆弧半径为3
L = 1.0;//设定单位长度

Point(1) = {L, -L, 0};//第四象限点1
Point(2) = {L, L, 0};//第一象限点2
Point(3) = {-L, L, 0};//第二象限点3
Point(4) = {-L, -L, 0};//第三象限点4
Point(5) = {-L + R, -L, 0};//圆弧起点5
Point(6) = {-L, -L + R, 0};//圆弧终点6
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};//圆弧中点7

Circle(1) = {5, 4, 7};//第一段圆弧，从点5到点7，圆心为点4
Circle(2) = {7, 4, 6};//第二段圆弧，从点7到点6，圆心为点4

Line(3) = {6, 3};//线3从6到3
Line(4) = {3, 2};//线4从3到2
Line(5) = {2, 1};//线5从2到1
Line(6) = {1, 5};//线6从1到5
Line(7) = {2, 7};//线7从2到7

Curve Loop(1) = {4, 7, 2, 3};
Plane Surface(1) = {1};

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;

Transfinite Surface{1};
Transfinite Surface{2};

Recombine Surface{1};
Recombine Surface{2};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;

// EOF

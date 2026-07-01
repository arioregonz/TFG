SetFactory("OpenCASCADE");

// --- Parámetros ---
pml = 0.05;    // grosor PML
h   = 0.01;    // tamaño elemento
Lx  = 1.00;   // dominio físico
Lz  = 1.00;
tx  = 0.35;   // inicio tronco en dominio físico
tz  = 0.35;
tw  = 0.30;   // tamaño tronco

// --- Número de elementos ---
n_pml  = 5;              // elementos en la capa PML (pml/h = 5)
n_tx1  = tx/h;           // elementos desde borde izq hasta tronco = 35
n_tx2  = tw/h;           // elementos en el tronco = 30
n_tx3  = (Lx-tx-tw)/h;  // elementos desde tronco hasta borde der = 35
n_tz1  = tz/h;           // igual en Z
n_tz2  = tw/h;
n_tz3  = (Lz-tz-tw)/h;

// =========================================================================
// PUNTOS — definición explícita de todos los puntos de la malla
// Coordenadas X: 0, pml, pml+tx, pml+tx+tw, pml+Lx, pml+Lx+pml
//              = 0, 0.05, 0.40,  0.70,       1.05,   1.10
// Coordenadas Z: 0, pml, pml+tz, pml+tz+tw, pml+Lz, pml+Lz+pml
//              = 0, 0.05, 0.40,  0.70,       1.05,   1.10
// =========================================================================

// X coords
x0=0; x1=pml; x2=pml+tx; x3=pml+tx+tw; x4=pml+Lx; x5=pml+Lx+pml;
// Z coords
z0=0; z1=pml; z2=pml+tz; z3=pml+tz+tw; z4=pml+Lz; z5=pml+Lz+pml;

// Fila z=z0 (6 puntos)
Point(1)  = {x0, z0, 0, h};
Point(2)  = {x1, z0, 0, h};
Point(3)  = {x2, z0, 0, h};
Point(4)  = {x3, z0, 0, h};
Point(5)  = {x4, z0, 0, h};
Point(6)  = {x5, z0, 0, h};

// Fila z=z1
Point(7)  = {x0, z1, 0, h};
Point(8)  = {x1, z1, 0, h};
Point(9)  = {x2, z1, 0, h};
Point(10) = {x3, z1, 0, h};
Point(11) = {x4, z1, 0, h};
Point(12) = {x5, z1, 0, h};

// Fila z=z2
Point(13) = {x0, z2, 0, h};
Point(14) = {x1, z2, 0, h};
Point(15) = {x2, z2, 0, h};
Point(16) = {x3, z2, 0, h};
Point(17) = {x4, z2, 0, h};
Point(18) = {x5, z2, 0, h};

// Fila z=z3
Point(19) = {x0, z3, 0, h};
Point(20) = {x1, z3, 0, h};
Point(21) = {x2, z3, 0, h};
Point(22) = {x3, z3, 0, h};
Point(23) = {x4, z3, 0, h};
Point(24) = {x5, z3, 0, h};

// Fila z=z4
Point(25) = {x0, z4, 0, h};
Point(26) = {x1, z4, 0, h};
Point(27) = {x2, z4, 0, h};
Point(28) = {x3, z4, 0, h};
Point(29) = {x4, z4, 0, h};
Point(30) = {x5, z4, 0, h};

// Fila z=z5
Point(31) = {x0, z5, 0, h};
Point(32) = {x1, z5, 0, h};
Point(33) = {x2, z5, 0, h};
Point(34) = {x3, z5, 0, h};
Point(35) = {x4, z5, 0, h};
Point(36) = {x5, z5, 0, h};

// =========================================================================
// LÍNEAS HORIZONTALES (por filas)
// =========================================================================
// Fila z0
Line(1)  = {1,2};   Line(2)  = {2,3};   Line(3)  = {3,4};
Line(4)  = {4,5};   Line(5)  = {5,6};
// Fila z1
Line(6)  = {7,8};   Line(7)  = {8,9};   Line(8)  = {9,10};
Line(9)  = {10,11}; Line(10) = {11,12};
// Fila z2
Line(11) = {13,14}; Line(12) = {14,15}; Line(13) = {15,16};
Line(14) = {16,17}; Line(15) = {17,18};
// Fila z3
Line(16) = {19,20}; Line(17) = {20,21}; Line(18) = {21,22};
Line(19) = {22,23}; Line(20) = {23,24};
// Fila z4
Line(21) = {25,26}; Line(22) = {26,27}; Line(23) = {27,28};
Line(24) = {28,29}; Line(25) = {29,30};
// Fila z5
Line(26) = {31,32}; Line(27) = {32,33}; Line(28) = {33,34};
Line(29) = {34,35}; Line(30) = {35,36};

// =========================================================================
// LÍNEAS VERTICALES (por columnas)
// =========================================================================
// Col x0
Line(31) = {1,7};   Line(32) = {7,13};  Line(33) = {13,19};
Line(34) = {19,25}; Line(35) = {25,31};
// Col x1
Line(36) = {2,8};   Line(37) = {8,14};  Line(38) = {14,20};
Line(39) = {20,26}; Line(40) = {26,32};
// Col x2
Line(41) = {3,9};   Line(42) = {9,15};  Line(43) = {15,21};
Line(44) = {21,27}; Line(45) = {27,33};
// Col x3
Line(46) = {4,10};  Line(47) = {10,16}; Line(48) = {16,22};
Line(49) = {22,28}; Line(50) = {28,34};
// Col x4
Line(51) = {5,11};  Line(52) = {11,17}; Line(53) = {17,23};
Line(54) = {23,29}; Line(55) = {29,35};
// Col x5
Line(56) = {6,12};  Line(57) = {12,18}; Line(58) = {18,24};
Line(59) = {24,30}; Line(60) = {30,36};

// =========================================================================
// SUPERFICIES (25 celdas = 5×5 grid)
// Nomenclatura: S(i,j) = celda en columna i, fila j (i,j de 1 a 5)
// =========================================================================
// Fila j=1 (z0 a z1) — PML inferior
Curve Loop(1)  = {1,36,-6,-31};   Plane Surface(1)  = {1};   // col1,fil1
Curve Loop(2)  = {2,41,-7,-36};   Plane Surface(2)  = {2};   // col2,fil1
Curve Loop(3)  = {3,46,-8,-41};   Plane Surface(3)  = {3};   // col3,fil1
Curve Loop(4)  = {4,51,-9,-46};   Plane Surface(4)  = {4};   // col4,fil1
Curve Loop(5)  = {5,56,-10,-51};  Plane Surface(5)  = {5};   // col5,fil1

// Fila j=2 (z1 a z2) — PML izq, dominio físico, PML der
Curve Loop(6)  = {6,37,-11,-32};  Plane Surface(6)  = {6};   // PML izq
Curve Loop(7)  = {7,42,-12,-37};  Plane Surface(7)  = {7};   // dominio
Curve Loop(8)  = {8,47,-13,-42};  Plane Surface(8)  = {8};   // dominio
Curve Loop(9)  = {9,52,-14,-47};  Plane Surface(9)  = {9};   // dominio
Curve Loop(10) = {10,57,-15,-52}; Plane Surface(10) = {10};  // PML der

// Fila j=3 (z2 a z3) — PML izq, aire, TRONCO, aire, PML der
Curve Loop(11) = {11,38,-16,-33}; Plane Surface(11) = {11};  // PML izq
Curve Loop(12) = {12,43,-17,-38}; Plane Surface(12) = {12};  // aire izq
Curve Loop(13) = {13,48,-18,-43}; Plane Surface(13) = {13};  // TRONCO
Curve Loop(14) = {14,53,-19,-48}; Plane Surface(14) = {14};  // aire der
Curve Loop(15) = {15,58,-20,-53}; Plane Surface(15) = {15};  // PML der

// Fila j=4 (z3 a z4) — PML izq, dominio físico, PML der
Curve Loop(16) = {16,39,-21,-34}; Plane Surface(16) = {16};  // PML izq
Curve Loop(17) = {17,44,-22,-39}; Plane Surface(17) = {17};  // dominio
Curve Loop(18) = {18,49,-23,-44}; Plane Surface(18) = {18};  // dominio
Curve Loop(19) = {19,54,-24,-49}; Plane Surface(19) = {19};  // dominio
Curve Loop(20) = {20,59,-25,-54}; Plane Surface(20) = {20};  // PML der

// Fila j=5 (z4 a z5) — PML superior
Curve Loop(21) = {21,40,-26,-35}; Plane Surface(21) = {21};  // col1,fil5
Curve Loop(22) = {22,45,-27,-40}; Plane Surface(22) = {22};  // col2,fil5
Curve Loop(23) = {23,50,-28,-45}; Plane Surface(23) = {23};  // col3,fil5
Curve Loop(24) = {24,55,-29,-50}; Plane Surface(24) = {24};  // col4,fil5
Curve Loop(25) = {25,60,-30,-55}; Plane Surface(25) = {25};  // col5,fil5

// =========================================================================
// TRANSFINITE — forzar elementos regulares en cada celda
// =========================================================================
// Líneas horizontales — número de nodos según la columna
// Col PML (x0-x1, x4-x5): n_pml+1 nodos
// Col dominio izq (x1-x2): n_tx1+1 nodos
// Col tronco (x2-x3): n_tx2+1 nodos
// Col dominio der (x3-x4): n_tx3+1 nodos

// Filas horizontales (todas iguales en cada fila)
Transfinite Curve {1,6,11,16,21,26}  = n_pml+1;   // x0→x1
Transfinite Curve {2,7,12,17,22,27}  = n_tx1+1;   // x1→x2
Transfinite Curve {3,8,13,18,23,28}  = n_tx2+1;   // x2→x3
Transfinite Curve {4,9,14,19,24,29}  = n_tx3+1;   // x3→x4
Transfinite Curve {5,10,15,20,25,30} = n_pml+1;   // x4→x5

// Líneas verticales
Transfinite Curve {31,32,33,34,35}   = n_pml+1;   // z0→z1 (PML inf)
Transfinite Curve {36,37,38,39,40}   = n_pml+1;   // misma col x1
Transfinite Curve {41,42,43,44,45}   = n_pml+1;
Transfinite Curve {46,47,48,49,50}   = n_pml+1;
Transfinite Curve {51,52,53,54,55}   = n_pml+1;
Transfinite Curve {56,57,58,59,60}   = n_pml+1;

// Corregir filas verticales por zona
// Fila j=1 (PML inf): ya puesto arriba como n_pml+1
// Fila j=2 (z1→z2): n_tz1+1
Transfinite Curve {32,37,42,47,52,57} = n_tz1+1;  // z1→z2
// Fila j=3 (z2→z3): n_tz2+1
Transfinite Curve {33,38,43,48,53,58} = n_tz2+1;  // z2→z3
// Fila j=4 (z3→z4): n_tz3+1
Transfinite Curve {34,39,44,49,54,59} = n_tz3+1;  // z3→z4
// Fila j=5 (PML sup): n_pml+1
Transfinite Curve {35,40,45,50,55,60} = n_pml+1;  // z4→z5

// Aplicar Transfinite a todas las superficies
Transfinite Surface {1,2,3,4,5,
                     6,7,8,9,10,
                     11,12,13,14,15,
                     16,17,18,19,20,
                     21,22,23,24,25};

// Recombinar todos los elementos en cuadriláteros
Recombine Surface {1,2,3,4,5,
                   6,7,8,9,10,
                   11,12,13,14,15,
                   16,17,18,19,20,
                   21,22,23,24,25};

// =========================================================================
// GRUPOS FÍSICOS
// M1 = aire + PML (todas las superficies excepto el tronco)
// M2 = tronco (superficie 13)
// Bordes exteriores del dominio total
// =========================================================================
Physical Surface("M1", 1) = {1,2,3,4,5,
                               6,7,8,9,10,
                               11,12,14,15,
                               16,17,18,19,20,
                               21,22,23,24,25};
Physical Surface("M2", 2) = {13};  // TRONCO

// Bordes exteriores (z=0, x=0, x=1.10, z=1.10)
Physical Line("Bottom", 1) = {1,2,3,4,5};     // z=0
Physical Line("Left",   2) = {31,32,33,34,35}; // x=0
Physical Line("Right",  3) = {56,57,58,59,60}; // x=1.10
Physical Line("Top",    4) = {26,27,28,29,30}; // z=1.10



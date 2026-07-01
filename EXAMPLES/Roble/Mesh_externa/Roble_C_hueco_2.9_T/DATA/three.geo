
SetFactory("OpenCASCADE");

// Parámetros
pml = 0.05;    // grosor PML
R   = 0.15;    // radio tronco (Ø30cm)
cx  = 0.55;    // centro X (dominio total 1.10m → centro en 0.55m)
cz  = 0.55;    // centro Z
Rh   = 0.0145;    // radio del hueco (Ø10cm)

// Dominio total 1.10 x 1.10 m (incluye PML)
Rectangle(1) = {0, 0, 0, 1.10, 1.10};

// Tronco cilíndrico
Disk(2) = {cx, cz, 0, R, R};
// Hueco cilíndrico
Disk(3) = {cx,cz, 0, Rh, Rh};

// Fragmentación — comparte nodos en la interfaz aire-tronco
BooleanFragments{ Surface{1}; Delete; }{ Surface{2,3}; Delete; }

// Control de calidad
Mesh.CharacteristicLengthMax = 0.01;
Mesh.CharacteristicLengthMin = 0.01;
Mesh.CharacteristicLengthMin = 0.01;
Mesh.RecombineAll             = 1;
Mesh.RecombinationAlgorithm   = 1;  // Blossom (mejor para cuadriláteros)
Mesh.Algorithm                = 6;  // Frontal-Delaunay (mejor para círculos)

// Grupos físicos
// IMPORTANTE: tras BooleanFragments los números cambian.
// Abrir en Gmsh GUI y verificar antes de mallar.
// Típicamente tras fragmentar Rectangle+Disk:
//   Surface{2} = tronco (disco)
//   Surface{3} = aire exterior (rectángulo con el disco recortado)
Physical Surface("M1", 1) = {3,4};   // aire — VERIFICAR EN GUI
Physical Surface("M2", 2) = {5};   // madera — VERIFICAR EN GUI

// Bordes absorbentes (exterior del dominio con PML)
Physical Line("Bottom", 1) = {1};  // VERIFICAR EN GUI
Physical Line("Left",   2) = {2};
Physical Line("Right",  3) = {3};
Physical Line("Top",    4) = {4};



SetFactory("OpenCASCADE");

// 1. Geometría (Dominio 1x1)
Rectangle(1) = {0, 0, 0, 1, 1};
Rectangle(2) = {0.35, 0.35, 0, 0.30, 0.30};

// 2. Fragmentación
BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete; }

// 3. CONTROL DE CALIDAD (Lo más importante)
// Esto dice: "haz los elementos de un tamaño máximo de 0.02"
// Al ser un valor pequeño, la malla será más fina automáticamente.
Mesh.CharacteristicLengthMax = 0.01; 

// Esto fuerza a Gmsh a crear cuadrángulos en lugar de triángulos, 
// que es lo que SPECFEM2D necesita para funcionar bien.
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 1;


// 4. Grupos Físicos (Verifica los números en la GUI de Gmsh)
// Tras BooleanFragments, mira en el árbol qué números tienen las superficies
Physical Surface("M1", 1) = {3}; 
Physical Surface("M2", 2) = {2};

// Bordes para las condiciones de contorno
Physical Line("Bottom", 1) = {1};
Physical Line("Left", 2) = {2};
Physical Line("Right", 3) = {3};
Physical Line("Top", 4) = {4};

// 5. Algoritmo de mallado
Mesh.Algorithm = 6;
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";


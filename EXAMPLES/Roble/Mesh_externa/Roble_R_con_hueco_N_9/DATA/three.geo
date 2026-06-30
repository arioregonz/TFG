SetFactory("OpenCASCADE"); 

// 1. Geometría (Dominio 1x1) 

// capa inferior
Rectangle(1) = {0.00,0.00,0,1.00,0.35};

// capa lateral izquierda
Rectangle(7) = {0.00,0.35,0,0.35,0.30};

// capa lateral derecha
Rectangle(8) = {0.65,0.35,0,0.35,0.30};

// capa superior
Rectangle(9) = {0.00,0.65,0,1.00,0.35};

// tronco inferior (parte inferior del 30x30)
 Rectangle(2) = {0.35,0.35,0,0.30,0.10};

// tronco superior
 Rectangle(3) = {0.35,0.55,0,0.30,0.10};

// tronco izquierdo
 Rectangle(4) = {0.35,0.45,0,0.10,0.10};

// tronco derecho
 Rectangle(5) = {0.55,0.45,0,0.10,0.10}; 

// tronco centro (anillo alrededor del hueco)
Rectangle(6) = {0.45,0.45,0,0.10,0.10};

// 2. Fragmentación
BooleanFragments{ Surface{1,7,8,9}; Delete; }{ Surface{2, 3, 4, 5, 6}; Delete; }

// 3. CONTROL DE CALIDAD (Lo más importante) 

// Esto dice: "haz los elementos de un tamaño máximo de 0.02" 
// Al ser un valor pequeño, la malla será más fina automáticamente. 

Mesh.CharacteristicLengthMax = 0.02; 

// Esto fuerza a Gmsh a crear cuadrángulos en lugar de triángulos, 
// que es lo que SPECFEM2D necesita para funcionar bien. 

Mesh.RecombineAll = 1; 
Mesh.RecombinationAlgorithm = 1; 

// 4. Grupos Físicos (Verifica los números en la GUI de Gmsh) 
// Tras BooleanFragments, mira en el árbol qué números tienen las superficies 

Physical Surface("M1",1) = {1,6,7,8,9}; 
Physical Surface("M2",2) = {2,3,4,5}; 

// Bordes para las condiciones de contorno 

Physical Line("Bottom", 1) = {1}; 
Physical Line("Left", 2) = {2}; 
Physical Line("Right", 3) = {3}; 
Physical Line("Top", 4) = {4}; 

// 5. Algoritmo de mallado

Mesh.Algorithm = 6; 

//+ Show "*"; 
//+ Show "*";

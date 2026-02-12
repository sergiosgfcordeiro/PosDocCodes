// GMSH script to generate a T3 structured mesh for a rectangular domain

// Geometry: Define rectangle dimensions
Lx = 48.0; // Length in x-direction
Ly1 = 44.0;  // Length in y-direction
Ly2 = 60.0;  // Length in y-direction

// Mesh resolution (number of divisions)
Nx = 65;  // Divisions along the x-axis
Ny = 65;   // Divisions along the y-axis

// Define points for the rectangle
Point(1) = {0, 0, 0, 1};   // Bottom-left
Point(2) = {Lx, Ly1, 0, 1};  // Bottom-right 
Point(3) = {Lx, Ly2, 0, 1}; // Top-left 
Point(4) = {0, Ly1, 0, 1};  // Top-right 

// Define the lines for the rectangle
Line(1) = {1, 2};  // Bottom side
Line(2) = {2, 3};  // Right side
Line(3) = {3, 4};  // Top side
Line(4) = {4, 1};  // Left side

// Define surfaces for the rectangles: cohesive interface and crack 
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define the mesh size
Mesh.CharacteristicLengthMax = 5.0;  // Max characteristic length for mesh elements

// Use the Transfinite command to create a structured mesh (grid)
Transfinite Line{1, 3} = Nx; // Horizontal lines (x-direction)
Transfinite Line{2, 4} = Ny; // Horizontal lines (x-direction)

// Set the transfinite mesh for the surface (structured elements)
Transfinite Surface {1} = {Nx, Ny};  // Apply structured mesh to the bottom panel

// Define physical surface
Physical Surface("cooksmembrane") = {1};  // Defines the panel

// Split the quadrilaterals into triangles
Mesh.Optimize = 1; // Enable optimization for better-quality mesh
Mesh.ElementOrder = 1;  // Set the element order for T3 (triangular elements)

// Set the mesh format
// Mesh.Format = 2;  // Mesh format 2.2
Mesh.MshFileVersion = 2.0;

// Specifies the output format 
// gmsh mesh.geo -format msh2 -o mesh.msh

// Generate the mesh
Mesh 2;

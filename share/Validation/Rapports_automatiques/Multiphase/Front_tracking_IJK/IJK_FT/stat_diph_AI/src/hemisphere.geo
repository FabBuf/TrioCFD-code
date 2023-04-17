// Gmsh project created on Mon May 13 19:51:40 2013

// Pas du maillage et rayon : 
// h=0.0001333333333333333;
h = 0.00045; // Je ne sais pas trop ce que h represente... Je ne crois pas que ce soit la longueur des arretes.
r = 0.004;

a = 1.0;
b = 1.0; 
c = 1.0;

// Equation : 
// (x/a)**2 + (y/b)**2 + (z/c)**2 = r**2
rx = r*a;
ry = r*b;
rz = r*c;

// Coordonnees du centre :
cx = 0.0001;
cy = 0.0075;
cz = 0.005;

// Points :
Point(1) = {cx   , cy   ,cz   , h}; // centre
Point(2) = {cx+rx, cy   ,cz   , h}; 
Point(3) = {cx   , cy+ry,cz   , h};
Point(4) = {cx   , cy   ,cz+rz, h};
Point(5) = {cx-rx*0.01, cy   ,cz   , h}; // On le laisse sur le plan pour dissimetrizer (cx-rx --> cx)
Point(6) = {cx   , cy-ry,cz   , h};
Point(7) = {cx   , cy   ,cz-rz, h};

// Point pour definir le grand axe de l'ellipse
Point(100) = {cx+r  , cy  ,cz  , h}; // Sur le grand axe Ox.
Point(200) = {cx  , cy+ry  ,cz  , h}; // Sur le grand axe Oy.

// Arc d'ellipse sur Ox : 
Ellipsis(1) = {2,1,100,3};
Ellipsis(2) = {3,1,100,5};
Ellipsis(3) = {5,1,100,6};
Ellipsis(4) = {6,1,100,2};
Ellipsis(5) = {2,1,100,7};
Ellipsis(6) = {7,1,100,5};
Ellipsis(7) = {5,1,100,4};
Ellipsis(8) = {4,1,100,2};
// A partir de la, elles sont dans le plan (yz) :
// donc le grand axe est forcement Oy :
// D'ou le point 200.
Ellipsis(9) = {6,1,200,7};
Ellipsis(10) = {7,1,200,3};
Ellipsis(11) = {3,1,200,4};
Ellipsis(12) ={4,1,200,6};

// Definition des contours par Line loop
Line Loop(1) = {1,11,8};
Line Loop(2) = {2,7,-11}; 
Line Loop(3) = {3,-12,-7}; 
Line Loop(4) = {4,-8,12}; 
Line Loop(5) = {5,10,-1}; 
Line Loop(6) = {-2,-10,6};
Line Loop(7) = {-3,-6,-9}; 
Line Loop(8) = {-4,9,-5}; 

// Ruled Surface : permet de d�finir des surfaces sph�riques � partir des contours. 
Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};
Ruled Surface(7) = {7};
Ruled Surface(8) = {8};

// Et Surface Loop : permet de d�finir la surface ferm�e engendr�e par les 8 surfaces d�finie ci-dessous. 
Surface Loop (1) = {1,2,3,4,5,6,7,8};

// Volume final : 
Volume (1) = {1};

// Quelques commandes : 
// gmsh hemisphere.geo & : Pour afficher la g�om�trie du r�sultat.
// gmsh hemisphere.geo -2 : Pour mailler la surface lat�rale.
// gmsh hemisphere.geo -3 : Pour mailler le volume de la sph�re.
// gmsh hemisphere.msh & : Pour afficher le maillage.

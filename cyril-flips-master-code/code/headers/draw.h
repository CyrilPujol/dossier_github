// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

#include "triangulation.h"
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>
typedef CGAL::Linear_cell_complex_for_combinatorial_map<2> LCC;
typedef LCC::Dart_handle Dart_handle_lcc;
typedef LCC::Point Point_lcc;

int plot(Triangulation tri, int draw_period_num=1);
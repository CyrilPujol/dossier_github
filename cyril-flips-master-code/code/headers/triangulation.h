// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

// Project headers
#ifndef TRIANGULATION_H
#define TRIANGULATION_H

// Include standard library utilities

#include <utility>
#include <vector>
#include <map>

// Include CGAL 

#include <CGAL/basic.h>
#include <CGAL/Combinatorial_map.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
//#include <CGAL/Arr_circle_segment_traits_2.h>
//#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_rational                          Rational;
typedef CGAL::Cartesian<Rational>                     Kernel;
typedef Kernel::Point_2                               Point;
typedef Kernel::Vector_2                              Vector ;   
typedef Kernel::Triangle_2                            Triangle ;   
typedef std::tuple<int,int>                           Isometrie;


/*!
\brief Class representing a triangulation of a hyperbolic surface.

The triangulation is represented by a combinatorial map. The geometry is represented by
a set of cross ratios : one for each edge of the triangulation. To construct a triangulation first construct a set of triangles by making calls to the method add_triangles.
Then glue the triangles together edge by edge by making calls to the method maked_edge. When calling make_edge you have to define the cross ratio of the edge so that the geometry
is defined throug this gluing process. In the end check the validity of the triangulation by calling the is_valid method. If the triangulation is to be drawn then it also contains the positions, in the
open unit disk, of the 3 vertices of one of its triangles. One of these 3 points is special in the sense that it is invariant under the flip operation.
*/

void print_status(std::string message);


class Info_vertex; //see far below

class Info_face; //see far below

class Triangulation{
private:
  /*!
  \brief Used by the combinatorial map of a triangulation to define the attribute types.

  Here we define only attributes for edges and they are cross ratios, .i.e complex numbers
  */


  struct Triangulation_cmap_item
  {
    template <class CMap>
    struct Dart_wrapper
    {
      typedef Isometrie   Dart_info; // the Isomerie of the dart
      typedef CGAL::Cell_attribute<CMap, Info_vertex * > Vertex_attrib; // Vertex attribute type 
      typedef CGAL::Cell_attribute<CMap, Info_face * > Face_attrib; // Face attribute type
      typedef std::tuple<Vertex_attrib,void,Face_attrib> Attributes; // Each dart will have attributes on faces and vertex
    };
  };

  typedef CGAL::Combinatorial_map<2, Triangulation_cmap_item>                         Combinatorial_map;
  typedef typename Combinatorial_map::template Attribute_handle<2>::type              Face_attribute_handle;
  typedef typename Combinatorial_map::template Attribute_handle<0>::type              Vertex_attribute_handle;


public:
  typedef typename Combinatorial_map::Dart                                            Dart;
  typedef typename Combinatorial_map::Dart_handle                                     DartHandle;
  typedef typename Combinatorial_map::Dart_range                                      DartRange;
  typedef typename Combinatorial_map::template One_dart_per_cell_range<0>             VertexRange;
  typedef typename Combinatorial_map::template One_dart_per_cell_range<1>             EdgeRange;
  typedef typename Combinatorial_map::template One_dart_per_cell_range<2>             FaceRange;

  // Initialisation

  Triangulation();

  // Creation

  void        clear        ();
  DartHandle  add_triangle (Info_face* info_f,Info_vertex* info_v1, Info_vertex* info_v2, Info_vertex* info_v3,Isometrie isom1,Isometrie isom2,Isometrie isom3 );
  void        make_edge    (DartHandle dart_1, DartHandle dart_2);
  void        set_period   (std::tuple<Rational,Rational> period_x, std::tuple<Rational,Rational> period_y);
  bool        is_edge_possible (DartHandle dart_1, DartHandle dart_2);
  bool        test_dart    (DartHandle dart);
  DartHandle  insert_point_in_triangle (DartHandle dart, Info_vertex* info_vi, Info_face* info_f1, Info_face* info_f2);

  // Access

  // -> darts manipulation
  DartHandle ccw(DartHandle dart);
  DartHandle cw(DartHandle dart);
  DartHandle opposite(DartHandle dart);
  

  int number_of_darts();

  DartRange   dart_range();
  VertexRange vertex_range();
  EdgeRange   edge_range();
  FaceRange   face_range();

  Info_face&   get_face_info  (DartHandle dart);
  Info_vertex& get_vertex_info(DartHandle dart);
  std::tuple<DartHandle,DartHandle,DartHandle> neighbors (DartHandle face);
  Rational     squared_length (DartHandle dart) ;
  Rational     pseudo_angle (Point p, Point q, Point r);
  std::tuple<Rational,Rational> angles_of_edge (DartHandle dh);
  Isometrie    get_isom       (DartHandle dart) ;
  bool         is_in_triangle (Info_vertex info_vi,Triangulation::DartHandle triangle);
  Point        real_position  (DartHandle dart) ;
  double       diameter       ();


  // Representation
  std::string isometrie_to_string (Isometrie isom); 
  void        print();

  // Isometries manipulation
  Isometrie isometrie_between  (Triangulation::DartHandle vertex1,Triangulation::DartHandle vertex2);
  Vector    isom_to_vect       (Isometrie isom);
  bool      is_more_left_bottom(Isometrie iso1, Isometrie iso2);
  Isometrie left_bottom_isom   (Isometrie iso1, Isometrie iso2, Isometrie iso3);
  Isometrie add_isom           (Isometrie iso1, Isometrie iso2 );
  Isometrie sub_isom           (Isometrie iso1, Isometrie iso2 );
  void      translate_triangle (DartHandle triangle, Isometrie offset);
  void      recenter();

  // Algorithmic
  void continuous_twist           (int direction=0,int coef=1);
  bool is_delaunay_point_triangle (DartHandle point, DartHandle triangle_1, DartHandle triangle_2, DartHandle triangle_3, Vector offset=Vector(0,0));
  bool is_delaunay();
  bool is_delaunay_edge           (DartHandle edge);
  void flip_edge                  (DartHandle edge_to_flip);

private:
  Combinatorial_map                 _cmap;
public:
  Vector period_x;
  Vector period_y;

};


 

class Info_vertex{
public:
  std::string label;
  Point position;
  std::list<Triangulation::DartHandle> touching_faces; 

  void init(std::string label, std::tuple<Rational,Rational> pos);

  std::string string_of_position();

  void add_touching_face(Triangulation::DartHandle dart);
};



class Info_face{
public:
  std::string label;

  void init(std::string label);
};

#endif
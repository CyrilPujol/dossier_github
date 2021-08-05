// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

/*
This file contains the implementation of the objects and functions declared in "triangulation.h".
*/

#include "../headers/triangulation.h"

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void print_status(std::string message){
    std::cout << message << std::endl;
}

//////////////////////////////////////////////////
/////////////// Initialisation ///////////////////
//////////////////////////////////////////////////

/*!
\brief Initializer of a vertex. Does not change the list of touching faces.
*/
void Info_vertex::init(std::string label, std::tuple<Rational,Rational> pos){
  Info_vertex::label = label;
  Info_vertex::position = Point(get<0>(pos), get<1>(pos)); 
}

/*!
\brief add the faces of dart to the list of touching faces.
*/
void Info_vertex::add_touching_face(Triangulation::DartHandle dart){
  Info_vertex::touching_faces.push_back(dart);
}

/*!
\brief return a string to print the position.
*/
std::string Info_vertex::string_of_position(){
  return "(" + std::to_string(CGAL::to_double(Info_vertex::position.x())) + "," + std::to_string(CGAL::to_double(Info_vertex::position.y())) + ")";
}

/*!
\brief Initializer of a face. Is now quite useless.
*/
void Info_face::init(std::string label){
  Info_face::label = label;
}

// Initialisation

/*!
\brief Initializer of the triangulation. Does not seem very usefull for now
*/
Triangulation::Triangulation(){
  clear();
}

//////////////////////////////////////////////////
// Creation //////////////////////////////////////
//////////////////////////////////////////////////

/*!
\brief Clears the triangulation.
*/
void Triangulation::clear(){
  _cmap.clear();
}

/*!
\brief Adds a triangle to the triangulation.
Constructs a lonely triangle and return one of its 3 darts. 
This method must be used when constructing the triangulation and may make the triangulation invalid.
*/
Triangulation::DartHandle Triangulation::add_triangle(Info_face* info_f,Info_vertex* info_v1, Info_vertex* info_v2, Info_vertex* info_v3,Isometrie isom1,Isometrie isom2,Isometrie isom3){
  Triangulation::DartHandle dh = _cmap.make_combinatorial_polygon(3);
  _cmap.set_attribute<2>(        dh, _cmap.create_attribute<2>( info_f ));
  _cmap.set_attribute<0>(        dh, _cmap.create_attribute<0>( info_v1));
  _cmap.set_attribute<0>(    cw(dh), _cmap.create_attribute<0>( info_v2));
  _cmap.set_attribute<0>(cw(cw(dh)), _cmap.create_attribute<0>( info_v3));
  Isometrie isom_left = left_bottom_isom(isom1,isom2,isom3);
  _cmap.info(dh) = sub_isom(isom1,isom_left);
  _cmap.info(cw(dh)) = sub_isom(isom2,isom_left);
  _cmap.info(ccw(dh)) = sub_isom(isom3,isom_left);
  return dh ;
}

/*!
\brief Makes an edge by gluing too darts.
This method must be used when constructing the triangulation.
*/
void Triangulation::make_edge(Triangulation::DartHandle dart_1, Triangulation::DartHandle dart_2){
  CGAL_assertion(_cmap.is_dart_used(dart_1) && _cmap.is_dart_used(dart_2) && _cmap.is_sewable<2>(dart_1, dart_2));
  _cmap.sew<2>(dart_1, dart_2);
}

/*!
\brief Set the period of the triangulation. Take tuples as input.
*/
void Triangulation::set_period(std::tuple<Rational,Rational> tuple_x, std::tuple<Rational,Rational> tuple_y){
  Point O = Point(0,0);
  Point x = Point(std::get<0>(tuple_x),std::get<1>(tuple_x));
  Point y = Point(std::get<0>(tuple_y),std::get<1>(tuple_y));
  period_x = Vector(O,x);
  period_y = Vector(O,y);
}

/*!
\brief Returns a bool.
*/
bool Triangulation::is_edge_possible(Triangulation::DartHandle dart_1, Triangulation::DartHandle dart_2){
  return ( _cmap.is_sewable<2>(dart_1, dart_2) );
}

/*!
\brief Used for DEBUG only
*/
bool Triangulation::test_dart(Triangulation::DartHandle dart_1){
  return (_cmap.is_dart_used(dart_1));
}

/*!
\brief insert a point corresponding to info_vi in the triangle of dart, i help name the faces
*/
Triangulation::DartHandle Triangulation::insert_point_in_triangle(Triangulation::DartHandle dart, Info_vertex* info_v, Info_face* info_f1, Info_face* info_f2){
  Isometrie iso1 = get_isom(    dart );
  Isometrie iso2 = get_isom( cw(dart));
  Isometrie iso3 = get_isom(ccw(dart));

  DartHandle newdart = _cmap.insert_cell_0_in_cell_2(dart);
  _cmap.set_attribute<0>(             newdart  , _cmap.create_attribute<0>( info_v));
  _cmap.set_attribute<0>(ccw(opposite(newdart)), _cmap.create_attribute<0>( info_v));
  _cmap.set_attribute<0>(opposite( cw(newdart)), _cmap.create_attribute<0>( info_v));

  //_cmap.set_attribute<2>(             newdart  , _cmap.create_attribute<2>( info_v));
  _cmap.set_attribute<2>(ccw(opposite(newdart)), _cmap.create_attribute<2>( info_f1));
  _cmap.set_attribute<2>(opposite( cw(newdart)), _cmap.create_attribute<2>( info_f2));

  _cmap.info(newdart)                = std::make_tuple(0,0);
  _cmap.info(ccw(opposite(newdart))) = std::make_tuple(0,0);
  _cmap.info(opposite(cw( newdart))) = std::make_tuple(0,0);

  _cmap.info(    opposite(newdart))  = iso1; 
  _cmap.info( opposite(ccw(opposite(newdart)))) = iso2;
  _cmap.info(cw(newdart)) = iso3;

  return newdart;
}

//////////////////////////////////////////////////
// Access ////////////////////////////////////////
//////////////////////////////////////////////////
/*!
\brief Returns the dart on counter-clockwise position from a given dart.
Equivalently Combinatorial_map::beta(dart, 1).
*/
Triangulation::DartHandle Triangulation::ccw(Triangulation::DartHandle dart){
  CGAL_assertion(!_cmap.template is_free<1>(dart));
  return _cmap.beta(dart, 1);
}

/*!
\brief Returns the dart on clockwise position from a given dart.
Equivalently Combinatorial_map::beta(dart, 0).
*/
Triangulation::DartHandle Triangulation::cw(Triangulation::DartHandle dart){
  CGAL_assertion(!_cmap.template is_free<0>(dart));
  return _cmap.beta(dart, 0);
}

/*!
\brief Returns the dart on opposite position from a given dart (along the same edge but not in the same triangle).
Equivalently Combinatorial_map::opposite(dart).
*/
Triangulation::DartHandle Triangulation::opposite(Triangulation::DartHandle dart){
  CGAL_assertion(!_cmap.template is_free<2>(dart));
  //return _cmap.beta(dart, 2);
  return _cmap.opposite(dart);
}

/*!
\brief Returns the number of darts of the triangulation.
*/
int Triangulation::number_of_darts(){
  return _cmap.number_of_darts();
}

/*!
\brief Returns a range over the darts of the triangulation.
*/
Triangulation::DartRange Triangulation::dart_range(){
  return _cmap.darts();
}

/*!
\brief Returns a range with handles to one dart per vertex.
*/
Triangulation::VertexRange Triangulation::vertex_range(){
  return _cmap.one_dart_per_cell<0>();
}

/*!
\brief Returns a range with handles to one dart per edge.
*/
Triangulation::EdgeRange Triangulation::edge_range(){
  return _cmap.one_dart_per_cell<1>();
}

/*!
\brief Returns a range with handles to one dart per face.
*/
Triangulation::FaceRange Triangulation::face_range(){
  return _cmap.one_dart_per_cell<2>();
}


/*!
\brief Returns a reference to the info associated to the face supported by the given dart
*/
Info_face& Triangulation::get_face_info(Triangulation::DartHandle dart){
  Face_attribute_handle attribute_handle = _cmap.attribute<2>(dart);
  CGAL_assertion(attribute_handle != 0);
  return * _cmap.info_of_attribute<2>(attribute_handle);
}
/*!
\brief Returns a reference to the info associated to the vertex supported by the given dart
*/
Info_vertex& Triangulation::get_vertex_info(Triangulation::DartHandle dart){
  Vertex_attribute_handle attribute_handle = _cmap.attribute<0>(dart);
  CGAL_assertion(attribute_handle != 0);
  return * _cmap.info_of_attribute<0>(attribute_handle);
}

/*!
/brief Returns a tuple containing the three neighboors faces of a triangle given by an edge. 
NB : The first face is the one next to the edge.
*/
std::tuple<Triangulation::DartHandle,Triangulation::DartHandle,Triangulation::DartHandle> Triangulation::neighbors(Triangulation::DartHandle face){
  Triangulation::DartHandle f1 = opposite(face);
  Triangulation::DartHandle f2 = opposite(cw(face));
  Triangulation::DartHandle f3 = opposite(ccw(face));
  return std::make_tuple(f1,f2,f3);
} 

/*!
\brief Get the length of the dart.
*/
Rational Triangulation::squared_length(Triangulation::DartHandle dart) {
  Point pos1 = real_position(    dart );
  Point pos2 = real_position(ccw(dart));
  Rational d = CGAL::squared_distance(pos1,pos2);
  return d;
}

/*!
\brief if angle is the (positive) angle pqr, returns f(angle) where f is monotonous  .
*/
Rational Triangulation::pseudo_angle (Point p, Point q, Point r){
  Vector v1 = p-q;
  Vector v2 = r-q;
  Rational dot_product = v1*v2;
  if (dot_product > 0)
  {return 1 - (dot_product * dot_product) / (v1.squared_length() * v2.squared_length());}
  else
  {return 1 + (dot_product * dot_product) / (v1.squared_length() * v2.squared_length());}
}

/*!
\brief returns the two PSEUDO angles opposite to the edge, in increasing order .
*/
std::tuple<Rational,Rational> Triangulation::angles_of_edge(Triangulation::DartHandle dh){
  DartHandle dh2 = cw(dh);
  DartHandle dh3 = cw(dh2);
  DartHandle dh_ = opposite(dh);
  DartHandle dh2_ = cw(dh_);
  DartHandle dh3_ = cw(dh2_);

  Rational angle1 = pseudo_angle(real_position(dh ),real_position(dh2 ),real_position(dh3 )); // la fonction CGAL n'existe pas
  Rational angle2 = pseudo_angle(real_position(dh_),real_position(dh2_),real_position(dh3_));  

  if (angle1 <= angle2)
  {return std::make_tuple(angle1,angle2);}
  else
  {return std::make_tuple(angle2,angle1);}
}


/*!
\brief Get the isometrie corresponding to the vertex pointed by the dart.
*/
Isometrie Triangulation::get_isom(Triangulation::DartHandle dart) {
  return   _cmap.info(dart)  ;
}

/*!
\brief Returns if the position of Info_vertex is in the triangle . 
*/
bool Triangulation::is_in_triangle(Info_vertex info_vi,Triangulation::DartHandle dh){
  auto point1 = real_position(    dh );
  auto point2 = real_position( cw(dh));
  auto point3 = real_position(ccw(dh));
  Triangle triangle(point1, point2, point3);
  auto p = info_vi.position;
  assert(triangle.bounded_side(p) !=0 );
  return triangle.bounded_side(p)>0;
}

Point Triangulation::real_position(DartHandle dh) {
  return get_vertex_info( dh ).position + isom_to_vect(get_isom( dh ));
}

double Triangulation::diameter(){
  double maxi = 0.;
  auto range1 = dart_range();
  for (auto it= range1.begin(); it!=range1.end(); ++it){
    auto range2 = dart_range();
    for (auto it2= range2.begin(); it2!=range2.end(); ++it2){
      maxi = std::max( CGAL::to_double(CGAL::squared_distance(real_position(it),real_position(it2))) , maxi );
    } 
  }
  return std::sqrt(maxi);
}


//////////////////////////////////////////////////
// Representation ////////////////////////////////
//////////////////////////////////////////////////

std::string Triangulation::isometrie_to_string(Isometrie isom){
  return  "(" + std::to_string(std::get<0>(isom)) +","+ std::to_string(std::get<1>(isom)) + ")" ;
}

/*!
\brief Prints the triangulation in the console.
*/
void Triangulation::print(){
  print_status("    _ P R I N T _");
  // 1. Make the combinatorial map display its characteristics (number of darts, etc...)
  _cmap.display_characteristics(std::cout);
  std::cout << std::endl;
  // 2. Print the verticies

  for (VertexRange::iterator it=vertex_range().begin(); it!=vertex_range().end(); ++it){
    Info_vertex& info = get_vertex_info(it);
    std::cout << info.label << " is in " << info.string_of_position() << std::endl;
  }

  // 3. Print the faces

  for (FaceRange::iterator it=face_range().begin(); it!=face_range().end(); ++it){
    Triangulation::DartHandle dh = it;
    std::cout << get_face_info(dh).label << " is made of vertex " 
        << get_vertex_info(    dh ).label << isometrie_to_string(get_isom(    dh )) <<", "
        << get_vertex_info( cw(dh)).label << isometrie_to_string(get_isom( cw(dh))) << ", " 
        << get_vertex_info(ccw(dh)).label << isometrie_to_string(get_isom(ccw(dh))) << std::endl;
  }

}
//////////////////////////////////////////////////
//////////// Isometries  /////////////////////////
//////////////////////////////////////////////////

/*!
\brief Returns the translation that goes from the first dart to the second. 
*/
Isometrie Triangulation::isometrie_between(Triangulation::DartHandle vertex1,Triangulation::DartHandle vertex2){
  Isometrie isom1 = get_isom(vertex1);
  Isometrie isom2 = get_isom(vertex2);
  return sub_isom(isom2,isom1) ;
}

/*!
\brief Transpose an isometrie as a vector. 
*/
Vector Triangulation::isom_to_vect(Isometrie isom){
  return (period_x * std::get<0>(isom)) + (period_y * std::get<1>(isom));
}

/*!
\brief returns if iso1 is left of (or right under) iso2. 
*/
bool Triangulation::is_more_left_bottom(Isometrie iso1, Isometrie iso2){
  return (std::get<0>(iso1) < std::get<0>(iso2) || ( std::get<0>(iso1) == std::get<0>(iso2) && std::get<1>(iso1) < std::get<1>(iso2) ) );
}

/*!
\brief Find the left-most isometrie which is the bottom-most one. 
*/
Isometrie Triangulation::left_bottom_isom(Isometrie iso1, Isometrie iso2, Isometrie iso3){
  if (is_more_left_bottom(iso1,iso2))
  {
    if (is_more_left_bottom(iso1,iso3))
    {
      return iso1;
    }
    else
    {
      return iso3;
    }
  }
  else
  {
    if (is_more_left_bottom(iso2,iso3))
    {
      return iso2;
    }
    else
    {
      return iso3;
    }
  } 
}


/*!
\brief Add isometries like vectors. 
*/
Isometrie Triangulation::add_isom(Isometrie iso1, Isometrie iso2 ){
  return std::make_tuple(std::get<0>(iso1)+std::get<0>(iso2),std::get<1>(iso1)+std::get<1>(iso2) );
}

/*!
\brief Substract isometries like vectors. 
*/
Isometrie Triangulation::sub_isom(Isometrie iso1, Isometrie iso2 ){
  return std::make_tuple(std::get<0>(iso1)-std::get<0>(iso2),std::get<1>(iso1)-std::get<1>(iso2) );
}

/*!
\brief Translate a triangle by a given offset. 
*/
void Triangulation::translate_triangle(Triangulation::DartHandle triangle,Isometrie offset){
  DartHandle a = triangle;  
  DartHandle b = cw(a);
  DartHandle c = cw(b);
  _cmap.info(a) = add_isom(_cmap.info(a),offset);
  _cmap.info(b) = add_isom(_cmap.info(b),offset);
  _cmap.info(c) = add_isom(_cmap.info(c),offset);
}

//////////////////////////////////////////////////
//// Algorithmic /////////////////////////////////
//////////////////////////////////////////////////

/*!
\brief apply a vertical/horizontal Dehn twist to the triangulation. 
*/
void Triangulation::continuous_twist(int direction, int coef){
  //direction is 0 for horizontal and 1 for vertical
  auto V = vertex_range();
  Isometrie iso;
  int layer;
  DartHandle dh;
  DartHandle dh_;
  for (auto it=V.begin(); it != V.end(); ++it){
    dh = ccw(cw(it)); // gives the const dart, thats not great ?
    dh_ = dh;
    iso = get_isom(dh);
    Point ancienne_position = get_vertex_info(dh).position;
    Point nouv_position;
    Rational ratio_x;
    Rational ratio_y;
    Isometrie decalage;
    auto range_dart = dart_range();
    switch (direction){
    case 0 : 
      ratio_x =  (Vector(Point(0,0),ancienne_position)* period_x)/period_x.squared_length() ;
      nouv_position = ancienne_position + coef*ratio_x*period_y;
      ratio_y =  (Vector(Point(0,0),nouv_position)* period_y)/period_y.squared_length() ;
      decalage = std::make_tuple(0, int(CGAL::to_double( ratio_y )) );
      do
      {
        Isometrie iso_ = get_isom(dh_);
        // On additionne la translation de base plus la translation induite par le point plus la translation induite par la translatino du point
        _cmap.info(dh_) = add_isom(add_isom(iso_,decalage),std::make_tuple(0,coef*std::get<0>(iso_)) ); 
        dh_= ccw(opposite(dh_));
      }
      while (! (dh_ == dh) );
      get_vertex_info(dh).position = nouv_position - (std::get<1>(decalage) * period_y);
      break;
    
    case 1 : 
      ratio_y =  (Vector(Point(0,0),ancienne_position)* period_y)/period_y.squared_length() ;
      nouv_position = ancienne_position + coef*ratio_y*period_x;
      ratio_x =  (Vector(Point(0,0),nouv_position)* period_x)/period_x.squared_length() ;
      decalage = std::make_tuple(int(CGAL::to_double( ratio_y )),0 );
      do
      {
        Isometrie iso_ = get_isom(dh_);
        // On additionne la translation de base plus la translation induite par le point plus la translation induite par la translatino du point
        _cmap.info(dh_) = add_isom(add_isom(iso_,decalage),std::make_tuple(coef*std::get<1>(iso_),0) ); 
        dh_= ccw(opposite(dh_));
      }
      while (! (dh_ == dh) );
      get_vertex_info(dh).position = nouv_position - (std::get<0>(decalage) * period_x);
      break;
    }
  }
  recenter();
}

/*!
\brief put every triangle in the right fondamental domain 
*/
void Triangulation::recenter(){
  auto F = face_range();
  for (auto it=F.begin(); it != F.end(); ++it){
    DartHandle dh1 = cw(it); // le cw s'assure de bien avoir un const dart
    DartHandle dh2 = cw(dh1);
    DartHandle dh3 = cw(dh2);
    Isometrie isoma = get_isom(dh1);
    Isometrie isomb = get_isom(dh2);
    Isometrie isomc = get_isom(dh3);
    Isometrie left_isom = left_bottom_isom(isoma,isomb,isomc);
    _cmap.info(dh1) = sub_isom(isoma,left_isom);
    _cmap.info(dh2) = sub_isom(isomb,left_isom);
    _cmap.info(dh3) = sub_isom(isomc,left_isom);
  }
}

/*!
\brief Returns if the point is ouside the interior of the bounded circle of the triangle. 
*/
bool Triangulation::is_delaunay_point_triangle(Triangulation::DartHandle point, Triangulation::DartHandle triangle_1,Triangulation::DartHandle triangle_2,Triangulation::DartHandle triangle_3, Vector offset){
  // triangle_1_2_3 are the vertex of a triangle, but this triangle is not already built, so one can not simply use a single dart.
  Point p = real_position(point) + offset;
  Point q = real_position(triangle_1);
  Point r = real_position(triangle_2);
  Point s = real_position(triangle_3);

  return( CGAL::side_of_bounded_circle(q,r,s,p) <= 0 ); //1 is inside, 0 is on, and -1 is outside the boundary
}

/*!
\brief Returns true whenever the edge supported by the dart edge_to_test is delaunay.

 edge->a __b        a__ b
       |\\ |   ->   | //|
      d|_\\|c      d|//_|c
*/
bool Triangulation::is_delaunay_edge(Triangulation::DartHandle edge){
  DartHandle a = edge;
  DartHandle b = cw(a);
  DartHandle c = opposite(a);
  DartHandle d = cw(c);
  Vector offset = isom_to_vect(get_isom(a))-isom_to_vect(get_isom(ccw(c)));
  bool check1 = is_delaunay_point_triangle(d,a,b,cw(b),offset);
  bool check2 = is_delaunay_point_triangle(b,c,d,cw(d),-offset);
  return ( check1 && check2 );
}

/*!
\brief Determines if the triangulation is delaunay (ie not flippable)
*/
bool Triangulation::is_delaunay(){
  //print_status("\nchecking if the triangulation is delaunay...");

  DartRange darts = dart_range();
  for (DartRange::iterator dart=darts.begin(); dart!=darts.end(); ++dart){
    if (! is_delaunay_edge(dart)){
      print_status("checked : the triangulation is not delaunay !!!");
      return false;
    }
  }
  //print_status("checked : the triangulation is delaunay");
  return true;
}


/*!
\brief Flips the edge supported by the dart edge_to_flip.
*/
void Triangulation::flip_edge(Triangulation::DartHandle edge_to_flip){
    //print_status("Going to flip an edge...");

    // First gather all the information needed
    //  d  b___c
    //   |\ \  |
    //   | \ \ | 
    //   |__\ \|
    //  f    e a
    //  a and d are the flipped darts


    DartHandle a = opposite(edge_to_flip);  //In case edge_to_flip is an iterator
    DartHandle b = ccw(a);
    DartHandle c = cw(a);

    DartHandle d = opposite(a);
    DartHandle e = ccw(d);
    DartHandle f = cw(d);

    Isometrie offset = isometrie_between(b,d);
    Isometrie offset_2 = isometrie_between(a,e);
    assert(offset==offset_2);

    translate_triangle(a,offset);

    Isometrie isoma = get_isom(a);
    Isometrie isomb = get_isom(b);
    Isometrie isomc = get_isom(c);
    Isometrie isomd = get_isom(d);
    Isometrie isome = get_isom(e);
    Isometrie isomf = get_isom(f);

    // Make the topological flip



    _cmap.unlink_beta<1>(a);
    _cmap.unlink_beta<1>(b);
    _cmap.unlink_beta<1>(c);

    _cmap.unlink_beta<1>(d);
    _cmap.unlink_beta<1>(e);
    _cmap.unlink_beta<1>(f);


    _cmap.link_beta<1>(b, a);
    _cmap.link_beta<0>(f, a); // This allows to conserve the information in f, not in a
    _cmap.link_beta<1>(f, b);

    _cmap.link_beta<1>(e, d);
    _cmap.link_beta<0>(c, d);// This allows to conserve the information in c, not in d
    _cmap.link_beta<1>(c, e);

    //actualize the isometries todo: deplacer avec un offset avant d'actualiser

    Isometrie left_isom = left_bottom_isom(isomc,isomb,isomf);
    _cmap.info(a) = sub_isom(isomc,left_isom);
    _cmap.info(b) = sub_isom(isomb,left_isom);
    _cmap.info(f) = sub_isom(isomf,left_isom);

    left_isom = left_bottom_isom(isomf,isome,isomc);
    _cmap.info(d) = sub_isom(isomf,left_isom) ;
    _cmap.info(e) = sub_isom(isome,left_isom) ;
    _cmap.info(c) = sub_isom(isomc,left_isom) ;

    Isometrie isoma_ = get_isom(a);
    Isometrie isomb_ = get_isom(b);
    Isometrie isomc_ = get_isom(c);
    Isometrie isomd_ = get_isom(d);
    Isometrie isome_ = get_isom(e);
    Isometrie isomf_ = get_isom(f);
    // CGAL_assertion(is_valid()); // TO BE REMOVED ONCE DEBUGGING IS FINISHED

    //print_status("...edge flipped");
}

// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

#include "../headers/triangulation.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*!
\brief get a printable string of the tuple (corresponding to a position or a translation)
*/
string tuple_to_string(tuple<Rational,Rational> t){
  return "(" + to_string(CGAL::to_double(get<0>(t))) + "," + to_string(CGAL::to_double(get<1>(t))) + ")" ;
}

/*!
\brief get the number of lines in the file
*/
int number_of_lines(std::string file_name) { 
    int number_lines = 0;
    std::string line;
    std::ifstream myfile(file_name);

    while (std::getline(myfile, line))
        ++number_lines;
    return number_lines;
}

/*!
\brief parse a line with two ints
*/
std::tuple<Rational,Rational> split_2_int(std::string s, std::string delimiter ){
  size_t pos = s.find(delimiter);
  Rational x = stod(s.substr(0, pos));
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  Rational y = stod(s.substr(0, pos));
  return make_tuple(x,y);
}
/*!
\brief parse a line with three ints
*/
std::tuple<int,int,int> split_3_int(std::string s, std::string delimiter ){
  size_t pos = s.find(delimiter);
  int x = stoi(s.substr(0, pos));
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int y = stoi(s.substr(0, pos));
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int z = stoi(s.substr(0, pos));
  return make_tuple(x,y,z);
}
/*!
\brief parse a line defining an isomertrie
*/
Isometrie split_isom(std::string s, std::string delimiter ){
  size_t pos = s.find(delimiter);
  std::string isom = s.substr(0, pos);
  if (isom != "isom") {std::cout<<"Fichier mal formate, les isometries doivent etre precedees de isom !\n";throw 1;}
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int x = stoi(s.substr(0, pos));
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int y = stoi(s.substr(0, pos));
  return make_tuple(x,y);
}
/*!
\brief parse a line defining the period
*/
std::tuple<Rational,Rational> split_period(std::string s, std::string delimiter ){
  size_t pos = s.find(delimiter);
  std::string period = s.substr(0, pos);
  if (period != "period") {std::cout<<"Fichier mal formate, les isometries doivent etre precedees de isom !\n";throw 1;}
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int x = stoi(s.substr(0, pos));
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  int y = stoi(s.substr(0, pos));
  return make_tuple(x,y);
}

/*!
\brief cut a line in two parts separated by the delimiter
*/
std::tuple<std::string,std::string> split_general(std::string s, std::string delimiter ){
  size_t pos = s.find(delimiter);
  std::string string1 = s.substr(0, pos);
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  std::string string2 = s.substr(0, pos);
  return make_tuple(string1,string2);
}

/*!
\brief initialize the settings
*/
std::map<std::string,std::string> create_map_settings(std::map<std::string,std::string> settings ,std::string file_name){
  std::string line;

  std::ifstream myfile (file_name);
  if (myfile.is_open())
  { 
    while ( getline (myfile,line) )
    {
      auto key_value = split_general(line,":");
      settings[std::get<0>(key_value)] = std::get<1>(key_value);
    }
    myfile.close();
  }
  else cout << "Unable to open file\n"; 
  return settings;
}

/*!
\brief initialize the verticies of the triangulation
*/
void initialise_vertices(Info_vertex *vertices, std::string file_name, int n) {
  std::string line;
  std::tuple<Rational,Rational> xy ;
  int i = 0;

  std::ifstream myfile (file_name);
  if (myfile.is_open())
  { 

    while ( getline (myfile,line) && i<n )
    {
      auto xy = split_2_int(line, " ");
      vertices[i].init("v"+to_string(i),xy);
      i++;
    }
    myfile.close();
  }
  else cout << "Unable to open file\n"; 
}

Triangulation initialise_period(Triangulation tri,std::string file_name) 
{
  std::ifstream myfile (file_name);
  if (myfile.is_open()){
    std::string line;
    // The file starts with two lines defining the period
    getline (myfile,line);
    std::tuple<Rational,Rational> period_x = split_period(line, " ");
    getline (myfile,line);
    std::tuple<Rational,Rational> period_y = split_period(line, " ");
    tri.set_period(period_x,period_y);
  }
  else{std::cout << "Unable to open file\n";}
  return tri;
}

/*!
\brief Initialize all the faces of the triangulation
*/
Triangulation initialise_faces(Info_face * faces, Info_vertex * vertices, Triangulation tri,int m,std::string file_name) {
  
  std::tuple<int,int,int> abc = make_tuple(0,0,0);
  int i = 0;

  std::ifstream myfile (file_name);
  if (myfile.is_open())
  {
    std::string line;
    while ( getline (myfile,line) )
    {
      // For each triangle, there are three isometries.
      abc = split_3_int(line, " ");
      getline (myfile,line);
      Isometrie isom1 = split_isom(line, " ");
      getline (myfile,line);
      Isometrie isom2 = split_isom(line, " ");
      getline (myfile,line);
      Isometrie isom3 = split_isom(line, " ");
      faces[i].init("f"+to_string(i));
      auto dh = tri.add_triangle(&faces[i],&vertices[get<0>(abc)],&vertices[get<1>(abc)],&vertices[get<2>(abc)],isom1,isom2,isom3);
      i++;
    }
    myfile.close();
  }
  else std::cout << "Unable to open file\n"; 
  return tri;
}

Triangulation initialise_triangulation_insertion(Info_vertex * vertices, Info_face * faces, int n, int m, Triangulation tri) {

  // On commence par ajouter ajouter un point en 0,0 et on initialise la triangulation avec lui.
  vertices[n-1].init("v*",std::make_tuple(0,0));
  for (int i = 0;i<m;i++){
    faces[i].init("f"+std::to_string(i));
  }
  auto dh0 = tri.add_triangle(&faces[m-1],&vertices[n-1],&vertices[n-1],&vertices[n-1],make_tuple(0,0),make_tuple(1,1),make_tuple(1,0));
  auto dh1 = tri.add_triangle(&faces[m-2],&vertices[n-1],&vertices[n-1],&vertices[n-1],make_tuple(0,0),make_tuple(0,1),make_tuple(1,1));
  tri.make_edge( tri.cw(dh0),        dh1 );
  tri.make_edge(tri.ccw(dh0),tri.cw( dh1));
  tri.make_edge(        dh0 ,tri.ccw(dh1));


  // On insère successivement chaque point de la triangulation. 
  // la triangulation initiale encecle le domaine fondamental, donc tout les points y seront inséré : pas de translation.
  Info_vertex info_vi;
  Triangulation::DartHandle dh;
  for (int i = 0; i < n-1; i++){ //le cas i=n-1 correspond à v*
    auto face_range = tri.face_range();
    for (auto it=face_range.begin(); it != face_range.end(); ++it){
      if (tri.is_in_triangle(vertices[i],it)){
        dh=tri.insert_point_in_triangle(it, &vertices[i], &faces[i+i], &faces[i+i+1]);
        // tout est bien collé en place ?
        break;
      }
    }
    //std::cout<<"just checking print "<<std::endl;
    //tri.print();
    
  }
  return tri;
}

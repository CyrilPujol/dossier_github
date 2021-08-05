// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

#include "triangulation.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

string tuple_to_string(tuple<Rational,Rational> t);

int number_of_lines(std::string file_name) ;

std::tuple<Rational,Rational> split_2_int(std::string s, std::string delimiter );
std::tuple<int,int,int> split_3_int(std::string s, std::string delimiter );
Isometrie split_isom(std::string s, std::string delimiter );
Isometrie split_period(std::string s, std::string delimiter );
std::tuple<std::string,std::string> split_general(std::string s, std::string delimiter );

std::map<std::string,std::string> create_map_settings(std::map<std::string,std::string> settings ,std::string file_name);

void initialise_vertices(Info_vertex *vertices,std::string file_name,int n) ;

Triangulation initialise_period(Triangulation tri,std::string file_name) ;

Triangulation initialise_faces(Info_face * faces, Info_vertex * vertices, Triangulation tri,int m,std::string file_name) ;

Triangulation initialise_triangulation_insertion(Info_vertex * vertices, Info_face * faces, int n, int m, Triangulation tri) ;


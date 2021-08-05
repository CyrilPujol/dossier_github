//#include "triangulation.cpp"
#include "../headers/draw.h"
#include "../headers/parsing.h"

using namespace std;

Triangulation tri;

/*!
\brief Compare two darts by length
*/
struct Cmp_length {
    bool operator()( Triangulation::DartHandle dh1,  Triangulation::DartHandle dh2) const
    {
        return (tri.squared_length(dh1) >= tri.squared_length(dh2)) ; // >= pour flipper les plus grands d'abbord
        // /!\ set operators are supposed to be strict, but here, two different darts can have the same length, 
        //     the "=" makes the set accept multiple darts of the same size, but it is kind of a cheat !
    }
};

struct Cmp_angle {
    bool operator()( Triangulation::DartHandle dh1,  Triangulation::DartHandle dh2) const
    {
        auto angles_1 = tri.angles_of_edge(dh1);
        auto angles_2 = tri.angles_of_edge(dh2);
        if (std::get<0>(angles_1) < std::get<0>(angles_2) || (std::get<0>(angles_1) == std::get<0>(angles_2) && std::get<1>(angles_1) <= std::get<1>(angles_2) ) )
        {return true;} // true+<= pour flipper les petits angles opposés, false et < sinon
        return false;
    }
};

/*!
\brief The waiting room orders the darts in a given manner
*/
//typedef std::stack<Triangulation::DartHandle>   waiting_room;
//typedef std::queue<Triangulation::DartHandle>   waiting_room;
typedef std::set<Triangulation::DartHandle,Cmp_angle>   waiting_room;

/*!
\brief Check if a dart is in the waiting room
*/
bool not_empty(waiting_room * room){
  return (! (*room).empty()); // for stack, queue and set
}

/*!
\brief Add a dart to the waiting room
*/
void add(waiting_room * room, Triangulation::DartHandle dh){
  //(*room).push(dh); // for stack and queue
  (*room).insert(dh); // for set 
}

/*!
\brief returns the next dart in the waiting room. It leaves the room.
*/
Triangulation::DartHandle next(waiting_room * room){
  Triangulation::DartHandle dh = * ((*room).begin()); // for set
  //Triangulation::DartHandle dh = (*room).front();// for queue 
  // (*room).top(); // for stack
  //
  //(*room).pop(); // for stack and queue
  (*room).erase((*room).begin()); // for set
  return(dh);
}

int main () {

  // This file contains all the settings, it should be correctly formatted, see README
  std::string file_settings = "data/settings.cin";

  

  std::map<std::string,std::string> settings = {};
  settings = create_map_settings(settings, file_settings);

  std::string file_points = settings["file_points"];
  std::string file_faces = settings["file_faces"];
  std::string file_period = settings["file_period"];
  std::string file_output;
  if (stoi(settings["output"])){ file_output = settings["file_output"]; }
  // récupérer les twists depuis un fichier externe si nécessaire
  std::ifstream myfile (settings["twists"]);
  if (myfile.is_open()) {  
    std::string line;
    
    getline (myfile,line);
    std::cout<<"number of twists : "<<line<<std::endl;
    settings["twists"] = line;
  }
  // récupérer la limite de points depuis un fichier externe si nécessaire
  std::ifstream myfile_ (settings["vertex_limit"]);
  if (myfile_.is_open()) {  // code moche, pour pouvoir récupérer les twists depuis un fichier externe
    std::string line;
    
    getline (myfile_,line);
    std::cout<<"vertex limit : "<<line<<std::endl;
    settings["vertex_limit"] = line;
  }

  tri = initialise_period(tri,file_period) ;

  int n = number_of_lines(file_points);
  if (stoi(settings["is_vertex_limit"])==1){n=std::min({n,stoi(settings["vertex_limit"])});}
  if (settings["input"]=="insertion"){n++;}
  Info_vertex vertices[n]; // set of all the vertices
  initialise_vertices(vertices,file_points,n);
  int m ;
  if (settings["input"]== "triangles")
  {m = (number_of_lines(file_faces))/4; /* 4 lignes par faces.*/ }
  else if (settings["input"]=="insertion")
  {m = 2*n; /* au début : 2 triangles puis chaque insertion remplace un triangle par trois */}
  Info_face faces[m]; // set of all information in triangles 

  if (settings["input"]=="triangles"){

    // create all the individual triangles
    tri = initialise_faces(faces,vertices,tri,m,file_faces);

    // add touching_faces, actually it's not very usefull
    auto T_raw = tri.face_range();
    for (Triangulation::FaceRange::iterator it = T_raw.begin(); !(it == T_raw.end()) ; ++it){
        Triangulation::DartHandle dh0 = it;
        Triangulation::DartHandle dh1 = tri.cw(dh0);
        Triangulation::DartHandle dh2 = tri.cw(dh1);
        tri.get_vertex_info(dh0).add_touching_face(dh0);
        tri.get_vertex_info(dh1).add_touching_face(dh1);
        tri.get_vertex_info(dh2).add_touching_face(dh2);
    }    

    tri.print();

    // sew the triangles together
    for (int i = 0; i<n; i++) // for every distinct vertex
    {
      Info_vertex vertex = vertices[i];
      auto touching_faces = vertex.touching_faces;
      for (auto it=touching_faces.begin(); it != touching_faces.end(); ++it) // for every face touching it
      {
        Triangulation::DartHandle dh0 = *it;
        Triangulation::DartHandle dh1 = tri.cw(dh0);
        Triangulation::DartHandle dh2 = tri.ccw(dh0);
        for (auto it2=touching_faces.begin(); it2 != touching_faces.end(); ++it2) // for every other face touching it
        {
          
          Triangulation::DartHandle dh0_ = *it2;
          if (tri.get_face_info(dh0).label < tri.get_face_info(dh0_).label)
          // face.label < face.label is an arbitrary order, it allows to compute once per pair of differents faces
          {   
            // The goal here is to check if the two triangles are supposed to be adjascent, and to glue them together if needed
            Isometrie delta = tri.isometrie_between(dh0,dh0_);
            Triangulation::DartHandle dh1_ = tri.cw(dh0_);
            Triangulation::DartHandle dh2_ = tri.ccw(dh0_);
            // This is very ugly and obfuscated it will probably be a problem later
            CGAL_assertion( tri.get_vertex_info(dh0).label == tri.get_vertex_info(dh0_).label );
            if (tri.get_vertex_info(dh2).label == tri.get_vertex_info(dh1_).label && tri.isometrie_between(dh2,dh1_) == delta)
            { tri.make_edge(dh0,dh1_);}
            if (tri.get_vertex_info(dh1).label == tri.get_vertex_info(dh2_).label && tri.isometrie_between(dh1,dh2_) == delta)
            { tri.make_edge(dh1,dh0_);}
            if (tri.get_vertex_info(dh1).label == tri.get_vertex_info(dh1_).label && tri.isometrie_between(dh1,dh1_) == delta)
            { tri.make_edge(dh1,dh1_);std::cout << "Warning : faces are not clockwise "<<tri.get_face_info(dh0).label<<tri.get_face_info(dh0_).label << std::endl;}
            if (tri.get_vertex_info(dh2).label == tri.get_vertex_info(dh2_).label && tri.isometrie_between(dh2,dh2_) == delta)
            { tri.make_edge(dh0,dh0_); std::cout << "Warning : faces are not clockwise "<<tri.get_face_info(dh0).label<<tri.get_face_info(dh0_).label << std::endl;}
          }
        }
      }
    }
    //tri.print();
  }
  
  else if (settings["input"]=="insertion"){

    tri = initialise_triangulation_insertion(vertices,faces,n,m, tri);

  }

  if ( stoi(settings["twists"])>0 && stoi(settings["draw"])>0 ){plot(tri,stoi(settings["draw_period_num"]));}
  
  // Twist the triangulation if needed
  if ( stoi(settings["twists"])>0 ){tri.continuous_twist(0,stoi(settings["twists"]));}

  if ( stoi(settings["draw"])==1 ){plot(tri,stoi(settings["draw_period_num"]));}

  double diameter;
  if (settings["input"]=="insertion"){diameter = std::sqrt( 1 + std::pow(1+stoi(settings["twists"]),2.0));}
  else{ diameter = tri.diameter();}

  auto E = tri.edge_range();
  //Now the flip part begins :
  waiting_room waiting_edges; // contains the edges to check and flip
  std::map<Triangulation::DartHandle, bool> is_waiting ; // say if a given dart is allready waiting
  int number_of_flips = 0;
  //add all the edges to the wating room
  for (auto it=E.begin(); it != E.end(); ++it) // for every edge 
  {
    add(&waiting_edges,it);
    is_waiting[it] = true ;
  }

  //While there are edges to check
  while (not_empty( &waiting_edges) ){
    Triangulation::DartHandle dh = next(&waiting_edges);
    is_waiting[dh] = false ;
    if (!tri.is_delaunay_edge(dh)){
      if ( stoi(settings["draw"])==2 ){plot(tri,stoi(settings["draw_period_num"]));}
      tri.flip_edge(dh);
      if (!tri.is_delaunay_edge(dh))
      {std::cout << "after a flip, the quadrilater should be Delaunay";throw 20;} 
      //actualiser les voisins
      auto neig1 = tri.neighbors(dh);
      auto neig2 = tri.neighbors(tri.opposite(dh));
      Triangulation::DartHandle edge1 = std::get<1>(neig1); //The four neighboors edges
      Triangulation::DartHandle edge2 = std::get<2>(neig1);
      Triangulation::DartHandle edge3 = std::get<1>(neig2);
      Triangulation::DartHandle edge4 = std::get<2>(neig2);
      // get<0> would give the edges we just swapped
      if (! is_waiting[edge1]) {add(&waiting_edges,edge1 ); is_waiting[edge1]=true;}; //Add them to the waiting room, if they are not allready here
      if (! is_waiting[edge2]) {add(&waiting_edges,edge2 ); is_waiting[edge2]=true;};
      if (! is_waiting[edge3]) {add(&waiting_edges,edge3 ); is_waiting[edge3]=true;};
      if (! is_waiting[edge4]) {add(&waiting_edges,edge4 ); is_waiting[edge4]=true;};

      //std::cout<< "edge flipped \n"; 
      number_of_flips ++;
      
    }
    else{
      //std::cout<< "edge not flipped \n"; 
    }
  }
  tri.is_delaunay() ;
  //tri.print();
  std::cout << "diameter :" << diameter << std::endl;
  std::cout << "number of edges flipped :" << number_of_flips << std::endl;

  if ( stoi(settings["draw"])>0 ){plot(tri,stoi(settings["draw_period_num"]));}

  if (stoi(settings["output"]))
  {
    std::ofstream out(file_output.c_str());
    if(out)    
    {out << diameter  <<" " << number_of_flips << std::endl;}
    else
    {cout << "ERREUR: Impossible d'ouvrir le fichier d'output." << endl;}
  }

  return 0; 
}
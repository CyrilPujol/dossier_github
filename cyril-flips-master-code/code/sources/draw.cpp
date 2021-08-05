// COPIRIGHT : Ce code est une propriété de l'INRIA
// AUTEUR : Cyril PUJOL

#include "../headers/draw.h"


int plot(Triangulation tri,int draw_period_num)
{
    LCC lcc;
    auto F = tri.face_range();

    for (auto it=F.begin(); it != F.end(); ++it) // for every face 
    {
        auto pos1 = tri.real_position(        it ) ;
        auto pos2 = tri.real_position( tri.cw(it));
        auto pos3 = tri.real_position(tri.ccw(it));

        double x1 = CGAL::to_double( pos1[0] );
        double y1 = CGAL::to_double( pos1[1] );
        double x2 = CGAL::to_double( pos2[0] );
        double y2 = CGAL::to_double( pos2[1] );
        double x3 = CGAL::to_double( pos3[0] );
        double y3 = CGAL::to_double( pos3[1] );

        for (int i=0; i<draw_period_num ;i++)
        {   

            for (int j=0; j<draw_period_num ;j++)
            {   
                Vector offset = (tri.period_x*i) + (tri.period_y*j);
                double offset_x = CGAL::to_double( offset[0] );
                double offset_y = CGAL::to_double( offset[1] );
                lcc.make_triangle(
                Point_lcc(x1 + offset_x ,y1 + offset_y ), 
                Point_lcc(x2 + offset_x ,y2 + offset_y ), 
                Point_lcc(x3 + offset_x ,y3 + offset_y ));
            }
        }
            
    }
    
    CGAL::draw(lcc,"Triangulation de Delaunay periodiques par flips");
    return EXIT_SUCCESS;
} 

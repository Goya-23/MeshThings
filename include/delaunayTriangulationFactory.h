#ifndef _DELAUNAY_TRIANGULATION_FACTORY_H_
#define _DELAUNAY_TRIANGULATION_FACTORY_H_

#include <memory>
#include "delaunayTriangulation.h"

class DelaunayTriangulationFactory
{
public:
    DelaunayTriangulationFactory() {}


    std::shared_ptr<DelaunayTriangulation> produce(char s);


private:
    std::shared_ptr<DelaunayTriangulation> m_dt;
};



#endif  //_DELAUNAY_TRIANGULATION_FACTORY_H_

#ifndef _DELAUNAY_TRIANGULATION_FACTORY_H_
#define _DELAUNAY_TRIANGULATION_FACTORY_H_

#include <memory>
#include <iostream>
#include "delaunayTriangulation.h"

class DelaunayTriangulationFactory
{
public:
    DelaunayTriangulationFactory() {}


    std::shared_ptr<DelaunayTriangulation> produce(char s);


};



#endif  //_DELAUNAY_TRIANGULATION_FACTORY_H_

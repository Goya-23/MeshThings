#include "delaunayTriangulationFactory.h"
#include "concreteDelaunayTriangulation.h"


std::shared_ptr<DelaunayTriangulation>
DelaunayTriangulationFactory::produce(char s)
{
    if(s == 'i' || s == 'I') {
        return std::make_shared<IncrementalDelaunayTriangulation>();
    } else if(s == 's' || s == 'S') {
        return std::make_shared<SweepLineDelaunayTriangulation>();
    } else {
        std::cout << "Invalid triangulation type" << std::endl;
        return nullptr;
    }
    
}
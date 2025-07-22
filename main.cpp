#include <iostream>
#include <memory>
#include "DelaunayTriangulationFactory.h"




int main(int argc, char** argv) {
    std::shared_ptr<PLCParser> plcParser = std::make_shared<PLCParser>();

    //i/I means incremental construction , s/S means sweep-line
    auto dtFactory = std::make_shared<DelaunayTriangulationFactory>();
    std::share_ptr<DelaunayTriangulation> dt = dtFactory->produce(argv[2][0]);
    dt->triangulate(plcParser->parse(argv[1]));


    return 0;
}
#ifndef _PLC_PARSER_H_
#define _PLC_PARSER_H_
#include <memory>
#include "plc.h"

class PLCParser
{
public:
    PLCParser() {}
    ~PLCParser() {}
    std::shared_ptr<PLC2D> parse(char* path);

};



#endif // _PLC_PARSER_H_ 
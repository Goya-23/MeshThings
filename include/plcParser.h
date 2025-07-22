#ifndef _PLC_PARSER_H_
#define _PLC_PARSER_H_
#include <memory>
#include "plc.h"

class PLCParser
{
public:
    PLCParser(char* path);
    ~PLCParser() {}
    std::shared_ptr<PLC2D> parse(char* path);
private:
    std::shared_ptr<PLC2D> plc;

}



#endif // _PLC_PARSER_H_ 
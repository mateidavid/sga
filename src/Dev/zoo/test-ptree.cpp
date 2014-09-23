#include <iostream>
#include "logger.hpp"
#include "ptree.hpp"

#define LOG_FACILITY "main"

using namespace std;


int main()
{
    log_l(info) << 42 << "\n";

    boost::property_tree::ptree t1;
    t1.put("t1_key", "val");
    log_l(info) << t1;

    ptree t2;
    t2.put("t2_key", "val");
    log_l(info) << t2;
}

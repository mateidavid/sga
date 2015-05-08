#include <cassert>
#include <iostream>
#include <tclap/CmdLine.h>
#include "logger.hpp"

#define TEST_FACILITY "test_facility"
#define ALT_FACILITY "alternate_facility"

using namespace std;

namespace global
{
    string program_name;
    string description =
        "Test Logger library. "
        "The program sends 5 log messages with decreasing priority (0=highest, 4=lowest) to 2 facilities. "
        "The level of '" TEST_FACILITY  "' is set explicitly, "
        "while the level of '" ALT_FACILITY "' is set via the default log level.";
    TCLAP::CmdLine cmd_parser(description);
    TCLAP::ValueArg< int > log_level        ("l", "log-level",         "Log level.",         false, 0, "int", cmd_parser);
    TCLAP::ValueArg< int > default_log_level("d", "default-log-level", "Deafult log level.", false, 0, "int", cmd_parser);
}

int real_main()
{
    cerr << "setting default log level to: " << global::default_log_level.getValue() << endl;
    Logger::set_default_level(global::default_log_level);

#define LOG_FACILITY TEST_FACILITY
    cerr << "setting " << LOG_FACILITY << " log level to: " << global::log_level.getValue() << endl;
    Logger::set_facility_level(LOG_FACILITY, global::log_level);
    int x = 0;
    LOG(error) << "this is an error message; x=" << x++ << endl;
    LOG(warning) << "this is a warning message; x=" << x++ << endl;
    LOG(info) << "this is an info message; x=" << x++ << endl;
    LOG(debug) << "this is a debug message; x=" << x++ << endl;
    LOG(4) << "this is a level 4 message; x=" << x++ << endl;
    cerr << "at the end, x=" << x << endl;
    assert(x == global::log_level + 1);
#undef LOG_FACILITY

#define LOG_FACILITY ALT_FACILITY
    x = 0;
    LOG(error) << "this is an error message; x=" << x++ << endl;
    LOG(warning) << "this is a warning message; x=" << x++ << endl;
    LOG(info) << "this is an info message; x=" << x++ << endl;
    LOG(debug) << "this is a debug message; x=" << x++ << endl;
    LOG(4) << "this is a level 4 message; x=" << x++ << endl;
    cerr << "at the end, x=" << x << endl;
    assert(x == global::default_log_level + 1);
#undef LOG_FACILITY

    cout << "ok\n";
    return 0;
}

int main(int argc, char* argv[])
{
    global::program_name = argv[0];
    global::cmd_parser.parse(argc, argv);
    real_main();
}

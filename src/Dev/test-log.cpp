#include <boost/program_options.hpp>
#include "Log.hpp"

#define LOG_FACILITY "test_facility"

using namespace std;


struct Program_Options
{
    int log_level;
    int default_log_level;
};

int real_main(const Program_Options& po)
{
    cerr << "setting default log level to: " << po.default_log_level << endl;
    detail::Log::default_level() = detail::Log::Level(po.default_log_level);

#define LOG_FACILITY "test_facility"
    cerr << "setting " << LOG_FACILITY << " log level to: " << po.log_level << endl;
    detail::Log::level(LOG_FACILITY) = detail::Log::Level(po.log_level);
    int x = 0;
    log_l(error) << "this is an error message; x=" << x++ << endl;
    log_l(warning) << "this is a warning message; x=" << x++ << endl;
    log_l(info) << "this is an info message; x=" << x++ << endl;
    log_l(debug) << "this is a debug message; x=" << x++ << endl;
    log_i(5) << "this is a level 5 message; x=" << x++ << endl;
    cerr << "at the end, x=" << x << endl;

#undef LOG_FACILITY
#define LOG_FACILITY "alternate_facility"
    x = 0;
    log_l(error) << "this is an error message; x=" << x++ << endl;
    log_l(warning) << "this is a warning message; x=" << x++ << endl;
    log_l(info) << "this is an info message; x=" << x++ << endl;
    log_l(debug) << "this is a debug message; x=" << x++ << endl;
    log_i(5) << "this is a level 5 message; x=" << x++ << endl;
    cerr << "at the end, x=" << x << endl;

    return 0;
}

int main(int argc, char* argv[])
{
    namespace bo = boost::program_options;
    Program_Options po;
    try
    {
        bo::options_description generic_opts_desc("Generic options");
        bo::options_description config_opts_desc("Configuration options");
        bo::options_description hidden_opts_desc("Hidden options");
        bo::options_description cmdline_opts_desc;
        bo::options_description visible_opts_desc("Allowed options");
        generic_opts_desc.add_options()
            ("help,h", "produce help message")
            ;
        config_opts_desc.add_options()
            ("log-level,l", bo::value<int>(&po.log_level)->default_value(0), "log level")
            ("default-log-level,d", bo::value<int>(&po.default_log_level)->default_value(0), "default log level")
            ;
        cmdline_opts_desc.add(generic_opts_desc).add(config_opts_desc).add(hidden_opts_desc);
        visible_opts_desc.add(generic_opts_desc).add(config_opts_desc);
        bo::variables_map vm;
        store(bo::command_line_parser(argc, argv).options(cmdline_opts_desc).run(), vm);
        notify(vm);
        if (vm.count("help"))
        {
            cout << visible_opts_desc;
            exit(EXIT_SUCCESS);
        }
    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
    real_main(po);
    return EXIT_SUCCESS;
}

#include <iostream>
#include "SGUtil.h"

using namespace std;


namespace global
{
    int verbosity;
    string prog_name;
}


void usage(ostream & os)
{
    os << "basic usage: " << global::prog_name
       << " -g <asqg_file>\n";
}


int main(int argc, char * argv[])
{
    global::prog_name = argv[0];

    string asqg_file;
    string out_file;

    char c;
    while ((c = getopt(argc, argv, "g:o:vh")) != -1)
    {
        istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
        case 'g':
            asqg_file = optarg;
            break;
        case 'o':
            out_file = optarg;
            break;
        case 'v':
            global::verbosity++;
            break;
        case 'h':
            usage(cout);
            exit(EXIT_SUCCESS);
        default:
            cerr << "unrecognized option: " << c << "\n";
            usage(cerr);
            exit(EXIT_FAILURE);
        }
    }
    if (optind != argc)
    {
        usage(cerr);
        exit(EXIT_FAILURE);
    }
    if (asqg_file.empty())
    {
        cerr << "asqg_file: not given\n";
        exit(EXIT_FAILURE);
    }
    if (out_file.empty())
    {
        cerr << "out_file: not given\n";
        exit(EXIT_FAILURE);
    }

    StringGraph* g_p = SGUtil::loadASQG(asqg_file, 0, true);
    clog << "loaded graph\n";

    g_p->writeASQG(out_file);

    delete g_p;

    return 0;
}

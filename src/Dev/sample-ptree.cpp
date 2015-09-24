#include <exception>
#include <set>
#include <vector>

#include "ptree.hpp"
#include <boost/property_tree/xml_parser.hpp>

using namespace std;

typedef boost::property_tree::ptree base_ptree;

void print_tree(ostream& os, const base_ptree& bpt)
{
    os << "XML:" << endl;
    boost::property_tree::write_xml(os, bpt);
    os << endl;
    os << "JSON:" << endl;
    try
    {
        boost::property_tree::write_json(os, bpt, false);
    }
    catch (exception& e)
    {
        os << "unprintable: " << e.what() << endl;
    }
}

#define PRINT_CMD(cmd) cout << "========== " BOOST_PP_STRINGIZE(cmd) << endl; cmd

int main()
{
    {
        base_ptree bpt;
        PRINT_CMD(bpt.put_value("val1"));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt;
        PRINT_CMD(bpt.put("", "val2"));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt;
        PRINT_CMD(bpt.put("key3", "val3"));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt;
        PRINT_CMD(bpt.put("key4_1.key4_2", "val4"));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt1;
        base_ptree bpt2;
        PRINT_CMD(bpt1.put("key1", "val1"));
        PRINT_CMD(bpt2.put("key2", "val2"));
        PRINT_CMD(bpt2.put_child("key3", bpt1));
        print_tree(cout, bpt2);
    }
    {
        base_ptree bpt;
        PRINT_CMD(auto& a = bpt.get_child(""));
        base_ptree bpt1;
        PRINT_CMD(bpt1.put_value("val1"));
        base_ptree bpt2;
        PRINT_CMD(bpt2.put("key2", "val2"));
        PRINT_CMD(a.push_back(std::make_pair("", bpt1)));
        PRINT_CMD(a.push_back(std::make_pair("", bpt2)));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt;
        PRINT_CMD(auto& a = bpt.get_child(""));
        base_ptree bpt1;
        PRINT_CMD(bpt1.put_value("val1"));
        base_ptree bpt2;
        PRINT_CMD(bpt2.put("key2", "val2"));
        PRINT_CMD(a.push_back(std::make_pair("bpt1", bpt1)));
        PRINT_CMD(a.push_back(std::make_pair("bpt2", bpt2)));
        print_tree(cout, bpt);
    }
    {
        base_ptree bpt;
        PRINT_CMD(auto& a = bpt.get_child(""));
        base_ptree bpt1;
        PRINT_CMD(bpt1.put_value("val1"));
        base_ptree bpt2;
        PRINT_CMD(bpt2.put("key2", "val2"));
        PRINT_CMD(a.push_back(std::make_pair("", bpt1)));
        PRINT_CMD(a.push_back(std::make_pair("", bpt2)));
        base_ptree bpt3;
        PRINT_CMD(bpt3.put_child("key0", bpt));
        print_tree(cout, bpt3);
    }
    {
        base_ptree bpt;
        PRINT_CMD(auto& a = bpt.get_child(""));
        base_ptree bpt1;
        PRINT_CMD(bpt1.put_value("val1"));
        base_ptree bpt2;
        PRINT_CMD(bpt2.put("key2", "val2"));
        PRINT_CMD(a.push_back(std::make_pair("bpt1", bpt1)));
        PRINT_CMD(a.push_back(std::make_pair("bpt2", bpt2)));
        base_ptree bpt3;
        PRINT_CMD(bpt3.put_child("key0", bpt));
        print_tree(cout, bpt3);
    }
    {
        ptree pt;
        set< int > s = { 1, 2, 3 };
        PRINT_CMD((pt.put("s1", range_to_ptree(s))));
        PRINT_CMD((pt.put("s2", s)));
        print_tree(cout, pt);
    }

}

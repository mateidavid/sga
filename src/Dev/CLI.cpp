#include "CLI.hpp"

namespace MAC
{

void cli(istream& is, ostream& os, const Graph&, const Hap_Map& hm)
{
    set< Mutation_Fact::index_type > unused_mut_set = Mutation_Fact::unused_set();
    set< Read_Chunk_Fact::index_type > unused_rc_set = Read_Chunk_Fact::unused_set();
    set< Read_Entry_Fact::index_type > unused_re_set = Read_Entry_Fact::unused_set();
    set< Contig_Entry_Fact::index_type > unused_ce_set = Contig_Entry_Fact::unused_set();
    set< Hap_Entry_Fact::index_type > unused_he_set = Contig_Entry_Fact::unused_set();
    set< Hap_Hop_Fact::index_type > unused_hh_set = Contig_Entry_Fact::unused_set();

    os << "Entering interactive mode. Type 'h' for help; 'q' or Ctrl-D to exit.\n";
    string line;
    while (getline(is, line))
    {
        istringstream iss(line + "\n");
        string cmd;
        if (not (iss >> cmd))
        {
            continue;
        }
        else if (cmd == "help" or cmd == "h")
        {
            os << "available commands:\n"
               << "  help|h : show this help\n"
               << "  quit|q : quit\n"
               << "  ptree|pt : print as ptree\n"
               << "  print|p : print object\n"
               << "  seq|s : print sequence\n"
               << "  outchunks|oc : print chunks leaving contig\n"
               << "  hops|hh : print haplotype hops\n"
                ;
        }
        else if (cmd == "quit" or cmd == "q")
        {
            break;
        }
        else if (cmd == "ptree" or cmd == "pt")
        {
            string tmp;
            size_t id;
            if (not (iss >> tmp >> id))
            {
                os << "use: ptree [(read|re)|(contig|ce)|(chunk|rc)] id\n";
            }
            else if (tmp == "read" or tmp == "re")
            {
                if (id >= Read_Entry_Fact::size() or unused_re_set.count(id) > 0)
                {
                    os << "Read_Entry[" << id << "]: not allocated\n";
                }
                else
                {
                    Read_Entry_CBPtr re_cbptr = Read_Entry_CBPtr::from_index(id);
                    os << "Read_Entry[" << id << "]: " << re_cbptr->to_ptree();
                }
            }
            else if (tmp == "contig" or tmp == "ce")
            {
                if (id >= Contig_Entry_Fact::size() or unused_ce_set.count(id) > 0)
                {
                    os << "Contig_Entry[" << id << "]: not allocated\n";
                }
                else
                {
                    Contig_Entry_CBPtr ce_cbptr = Contig_Entry_CBPtr::from_index(id);
                    os << "Contig_Entry[" << id << "]: " << ce_cbptr->to_ptree();
                }
            }
            else if (tmp == "chunk" or tmp == "rc")
            {
                if (id >= Read_Chunk_Fact::size() or unused_rc_set.count(id) > 0)
                {
                    os << "Read_Chunk[" << id << "]: not allocated\n";
                }
                else
                {
                    Read_Chunk_CBPtr rc_cbptr = Read_Chunk_CBPtr::from_index(id);
                    os << "Read_Chunk[" << id << "]: " << rc_cbptr->to_ptree();
                }
            }
            else
            {
                os << "invalid object to print: " << tmp << "\n";
            }
        }
        else if (cmd == "print" or cmd == "p")
        {
            string tmp;
            size_t id;
            if (not (iss >> tmp >> id))
            {
                os << "use: prettyprint [(read|re)|(contig|ce)|(chunk|rc)|(hap|he)|(hop|hh)] id\n";
            }
            else if (tmp == "read" or tmp == "re")
            {
                if (id >= Read_Entry_Fact::size() or unused_re_set.count(id) > 0)
                {
                    os << "Read_Entry:" << id << ": not allocated\n";
                }
                else
                {
                    Read_Entry_CBPtr re_cbptr = Read_Entry_CBPtr::from_index(id);
                    os << "Read_Entry: " << id
                       << " name: " << re_cbptr->name()
                       << " start: " << re_cbptr->start()
                       << " len: " << re_cbptr->end() - re_cbptr->start() << "\n";
                    for (auto rc_cbptr : re_cbptr->chunk_cont() | referenced)
                    {
                        os << "  " << Read_Chunk::to_string(rc_cbptr, true, true) << "\n";
                    }
                }
            }
            else if (tmp == "contig" or tmp == "ce")
            {
                if (id >= Contig_Entry_Fact::size() or unused_ce_set.count(id) > 0)
                {
                    os << "Contig_Entry:" << id << ": not allocated\n";
                }
                else
                {
                    Contig_Entry_CBPtr ce_cbptr = Contig_Entry_CBPtr::from_index(id);
                    os << "Contig_Entry: " << id
                       << " len: " << ce_cbptr->len() << "\n";
                    for (auto rc_cbptr : ce_cbptr->chunk_cont() | referenced)
                    {
                        os << "  " << Read_Chunk::to_string(rc_cbptr, false, true) << "\n";
                    }
                    for (auto mut_cbptr : ce_cbptr->mut_cont() | referenced)
                    {
                        os << "  ";
                        Mutation::to_stream(os, mut_cbptr, ce_cbptr);
                        os << endl;
                    }
                }
            }
            else if (tmp == "chunk" or tmp == "rc")
            {
                if (id >= Read_Chunk_Fact::size() or unused_rc_set.count(id) > 0)
                {
                    os << "Read_Chunk:" << id << ": not allocated\n";
                }
                else
                {
                    Read_Chunk_CBPtr rc_cbptr = Read_Chunk_CBPtr::from_index(id);
                    os << "  " << Read_Chunk::to_string(rc_cbptr) << "\n";
                }
            }
            else if (tmp == "hap" or tmp == "he")
            {
                if (id >= Hap_Entry_Fact::size() or unused_he_set.count(id) > 0)
                {
                    os << "Hap_Entry:" << id << ": not allocated\n";
                }
                else
                {
                    Hap_Entry_CBPtr he_cbptr = Hap_Entry_CBPtr::from_index(id);
                    os << "Hap_Entry: " << id << endl;
                    for (auto hh_cbptr : he_cbptr->hh_cont() | referenced)
                    {
                        os << "  ";
                        Hap_Hop::to_stream(os, hh_cbptr);
                        os << endl;
                    }
                }
            }
            else if (tmp == "hop" or tmp == "hh")
            {
                if (id >= Hap_Hop_Fact::size() or unused_hh_set.count(id) > 0)
                {
                    os << "Hap_Hop:" << id << ": not allocated\n";
                }
                else
                {
                    Hap_Hop_CBPtr hh_cbptr = Hap_Hop_CBPtr::from_index(id);
                    os << "  ";
                    Hap_Hop::to_stream(os, hh_cbptr);
                    os << endl;
                }
            }
            else
            {
                os << "invalid object to print: " << tmp << "\n";
            }
        }
        else if (cmd == "seq" or cmd == "s")
        {
            string tmp;
            size_t id;
            if (not (iss >> tmp >> id))
            {
                os << "use: seq [(read|re)|(contig|ce)|(chunk|rc)] id\n";
            }
            else if (tmp == "read" or tmp == "re")
            {
                if (id >= Read_Entry_Fact::size() or unused_re_set.count(id) > 0)
                {
                    os << "Read_Entry:" << id << ": not allocated\n";
                }
                else
                {
                    Read_Entry_CBPtr re_cbptr = Read_Entry_CBPtr::from_index(id);
                    os << "Read_Entry: " << id
                       << " name: " << re_cbptr->name()
                       << " start: " << re_cbptr->start()
                       << " len: " << re_cbptr->end() - re_cbptr->start() << "\n"
                       << re_cbptr->get_seq() << "\n";
                }
            }
            else if (tmp == "contig" or tmp == "ce")
            {
                if (id >= Contig_Entry_Fact::size() or unused_ce_set.count(id) > 0)
                {
                    os << "Contig_Entry:" << id << ": not allocated\n";
                }
                else
                {
                    Contig_Entry_CBPtr ce_cbptr = Contig_Entry_CBPtr::from_index(id);
                    os << "Contig_Entry: " << id
                       << " len: " << ce_cbptr->len() << endl
                       << ce_cbptr->seq() << endl;
                }
            }
            else if (tmp == "chunk" or tmp == "rc")
            {
                if (id >= Read_Chunk_Fact::size() or unused_rc_set.count(id) > 0)
                {
                    os << "Read_Chunk:" << id << ": not allocated\n";
                }
                else
                {
                    Read_Chunk_CBPtr rc_cbptr = Read_Chunk_CBPtr::from_index(id);
                    os << rc_cbptr->get_seq() << "\n";
                }
            }
            else
            {
                os << "invalid sequence to print: " << tmp << "\n";
            }
        }
        else if (cmd == "outchunks" or cmd == "oc")
        {
            size_t id;
            int c_right;
            int unmappable_policy = 1;
            int ignore_threshold = 0;
            iss >> id >> c_right;
            if (iss.eof())
            {
                os << "use: outchunks ce c_right [unmappable_policy [ignore_thres]]\n";
                continue;
            }
            iss >> unmappable_policy >> ignore_threshold; // might be absent
            if (id >= Contig_Entry_Fact::size() or unused_ce_set.count(id) > 0)
            {
                os << "Contig_Entry:" << id << ": not allocated\n";
                continue;
            }
            Contig_Entry_CBPtr ce_cbptr = Contig_Entry_CBPtr::from_index(id);
            auto cks = ce_cbptr->out_chunks_dir(c_right, unmappable_policy, ignore_threshold);
            for (auto t : cks)
            {
                Contig_Entry_CBPtr ce_next_cbptr;
                bool same_orientation;
                set< Read_Chunk_CBPtr > rc_cbptr_cont;
                tie(ce_next_cbptr, same_orientation) = move(t.first);
                rc_cbptr_cont = move(t.second);
                os << "  ce " << setw(9) << left << ce_next_cbptr.to_int()
                   << " same_orientation " << same_orientation
                   << " rc";
                for (const auto& rc_cbptr : rc_cbptr_cont)
                {
                    os << " " << rc_cbptr.to_int();
                }
                os << "\n";
            }
        }
        else if (cmd == "hops" or cmd == "hh")
        {
            size_t id;
            int c_right;
            bool is_endpoint;
            iss >> id;
            if (iss.eof())
            {
                os << "use: hops [ ce c_right | mut ]" << endl;
                continue;
            }
            iss >> c_right;
            is_endpoint = not iss.eof();
            Allele_Anchor anchor;
            if (is_endpoint)
            {
                if (id >= Contig_Entry_Fact::size() or unused_ce_set.count(id) > 0)
                {
                    os << "Contig_Entry:" << id << ": not allocated\n";
                    continue;
                }
                Contig_Entry_CBPtr ce_cbptr = Contig_Entry_CBPtr::from_index(id);
                anchor = Allele_Anchor(ce_cbptr, c_right);
            }
            else // is_mutation
            {
                if (id >= Mutation_Fact::size() or unused_mut_set.count(id) > 0)
                {
                    os << "Mutation:" << id << ": not allocated\n";
                    continue;
                }
                Mutation_CBPtr mut_cbptr = Mutation_CBPtr::from_index(id);
                anchor = Allele_Anchor(mut_cbptr);
            }
            auto p = hm.hh_set().equal_range(anchor);
            for (auto it = p.first; it != p.second; ++it)
            {
                os << "  ";
                Hap_Hop::to_stream(os, &*it);
                os << endl;
            }
        }
        else
        {
            os << "unknown command: " << cmd << "\n";
        }
    }
}

} // namespace MAC

#include "Mutation.hpp"

using namespace std;
using namespace MAC;


int main()
{
    Mutation_Fact mut_fact;
    Mutation_Chunk_Adapter_Fact mca_fact;

    // ctor: default
    Mutation_BPtr m1_ptr = Mutation_Fact::new_elem();
    assert(m1_ptr->rf_start() == 0);
    assert(m1_ptr->rf_end() == 0);
    assert(m1_ptr->rf_len() == 0);
    assert(m1_ptr->seq_len() == 0);
    assert(m1_ptr->have_seq());
    assert(m1_ptr->seq().size() == 0);
    assert(not m1_ptr->is_snp());
    assert(not m1_ptr->is_ins());
    assert(not m1_ptr->is_del());
    assert(m1_ptr->is_empty());

    // ctor: no alternate seq or seq_len
    Mutation_BPtr m2_ptr = Mutation_Fact::new_elem(3, 1);
    assert(m2_ptr->rf_start() == 3);
    assert(m2_ptr->rf_end() == 4);
    assert(m2_ptr->rf_len() == 1);
    assert(m2_ptr->seq_len() == 0);
    assert(m2_ptr->have_seq());
    assert(m2_ptr->seq().size() == 0);
    assert(not m2_ptr->is_snp());
    assert(not m2_ptr->is_ins());
    assert(m2_ptr->is_del());
    assert(not m2_ptr->is_empty());

    // ctor: alternate sequence
    Mutation_BPtr m3_ptr = Mutation_Fact::new_elem(5, 0, string("C"));
    assert(m3_ptr->rf_start() == 5);
    assert(m3_ptr->rf_end() == 5);
    assert(m3_ptr->rf_len() == 0);
    assert(m3_ptr->seq_len() == 1);
    assert(m3_ptr->have_seq());
    assert(m3_ptr->seq() == string("C"));
    assert(not m3_ptr->is_snp());
    assert(m3_ptr->is_ins());
    assert(not m3_ptr->is_del());
    assert(not m3_ptr->is_empty());

    // ctor:: alternate seq_len only
    Mutation_BPtr m4_ptr = Mutation_Fact::new_elem(8, 1, 1);
    assert(m4_ptr->rf_start() == 8);
    assert(m4_ptr->rf_end() == 9);
    assert(m4_ptr->rf_len() == 1);
    assert(m4_ptr->seq_len() == 1);
    assert(not m4_ptr->have_seq());
    assert(m4_ptr->is_snp());
    assert(not m4_ptr->is_ins());
    assert(not m4_ptr->is_del());
    assert(not m4_ptr->is_empty());

    cout << "m1:" << m1_ptr->to_ptree();
    cout << "m2:" << m2_ptr->to_ptree();
    cout << "m3:" << m3_ptr->to_ptree();
    cout << "m4:" << m4_ptr->to_ptree();

    Mutation_Fact::del_elem(m1_ptr);
    Mutation_Fact::del_elem(m2_ptr);
    Mutation_Fact::del_elem(m3_ptr);
    Mutation_Fact::del_elem(m4_ptr);

    return EXIT_SUCCESS;
}

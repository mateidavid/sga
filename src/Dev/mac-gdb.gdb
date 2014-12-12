source ~/code/interval-tree/include/boost/intrusive/itree.gdb
source ~/code/sga/src/Dev/factory.gdb

python

def get_first_subclass(v):
    res = ''
    if hasattr(v, 'qualifiers') and v.qualifiers:
        res += '(' + v.qualifiers + ')'
    res += str(v.cast(v.type.fields()[0].type))
    return res

def identity(v):
    return v

boost_print.inner_type[("MAC::detail::Mutation_ITree_Value_Traits", "node_traits")] = "MAC::detail::Mutation_ITree_Node_Traits"
boost_print.static_method[("MAC::detail::Mutation_ITree_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost_print.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost_print.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost_print.add_trivial_printer('MAC::Mutation_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Mutation_Ptr_List_Value_Traits", "node_traits")] = "MAC::detail::Mutation_Ptr_List_Node_Traits"
boost_print.static_method[("MAC::detail::Mutation_Ptr_List_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Mutation_Ptr_List_Node_Traits", "get_next")] = lambda v: v['_mut_ptr_next']
boost_print.add_trivial_printer('MAC::Mutation_Ptr_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Read_Chunk_ITree_Value_Traits", "node_traits")] = "MAC::detail::Read_Chunk_ITree_Node_Traits"
boost_print.static_method[("MAC::detail::Read_Chunk_ITree_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_left")] = lambda v: v['_ce_l_child']
boost_print.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_right")] = lambda v: v['_ce_r_child']
boost_print.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_parent")] = lambda v: v['_ce_parent']
boost_print.add_trivial_printer('MAC::Read_Chunk_CE_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Read_Chunk_Set_Value_Traits", "node_traits")] = "MAC::detail::Read_Chunk_Set_Node_Traits"
boost_print.static_method[("MAC::detail::Read_Chunk_Set_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_left")] = lambda v: v['_re_l_child']
boost_print.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_right")] = lambda v: v['_re_r_child']
boost_print.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_parent")] = lambda v: v['_re_parent']
boost_print.add_trivial_printer('MAC::Read_Chunk_RE_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Read_Chunk_Ptr_List_Value_Traits", "node_traits")] = "MAC::detail::Read_Chunk_Ptr_List_Node_Traits"
boost_print.static_method[("MAC::detail::Read_Chunk_Ptr_List_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Read_Chunk_Ptr_List_Node_Traits", "get_next")] = lambda v: v['_chunk_ptr_next']
boost_print.add_trivial_printer('MAC::Read_Chunk_Ptr_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Read_Entry_Set_Value_Traits", "node_traits")] = "MAC::detail::Read_Entry_Set_Node_Traits"
boost_print.static_method[("MAC::detail::Read_Entry_Set_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost_print.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost_print.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost_print.add_trivial_printer('MAC::Read_Entry_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Contig_Entry_List_Value_Traits", "node_traits")] = "MAC::detail::Contig_Entry_List_Node_Traits"
boost_print.static_method[("MAC::detail::Contig_Entry_List_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Contig_Entry_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost_print.add_trivial_printer('MAC::Contig_Entry_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Hap_Entry_List_Value_Traits", "node_traits")] = "MAC::detail::Hap_Entry_List_Node_Traits"
boost_print.static_method[("MAC::detail::Hap_Entry_List_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Hap_Entry_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost_print.add_trivial_printer('MAC::Hap_Entry_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Hap_Hop_List_Value_Traits", "node_traits")] = "MAC::detail::Hap_Hop_List_Node_Traits"
boost_print.static_method[("MAC::detail::Hap_Hop_List_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Hap_Hop_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost_print.add_trivial_printer('MAC::Hap_Hop_Cont', get_first_subclass)

boost_print.inner_type[("MAC::detail::Hap_Hop_Set_Value_Traits", "node_traits")] = "MAC::detail::Hap_Hop_Set_Node_Traits"
boost_print.static_method[("MAC::detail::Hap_Hop_Set_Value_Traits", "to_value_ptr")] = identity
boost_print.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost_print.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost_print.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost_print.add_trivial_printer('MAC::Hap_Hop_Set', get_first_subclass)

end

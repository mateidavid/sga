#!/usr/bin/env python
source itree.gdb
source factory.gdb
python

def get_first_subclass(v):
    qual = ''
    if hasattr(v, 'qualifiers'):
        qual = v.qualifiers
    tv = v.cast(v.type.fields()[0].type)
    if 'c' in qual:
        tv = tv.cast(tv.type.const())
    if 'v' in qual:
        tv = tv.cast(tv.type.volatile())
    if '&' in qual:
        tv = tv.cast(tv.type.reference())
    return tv

def identity(v):
    return v

boost.add_trivial_printer('rc_sequence::Sequence', get_first_subclass)

boost.static_method[("MAC::detail::Mutation_ITree_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost.static_method[("MAC::detail::Mutation_ITree_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost.add_trivial_printer('MAC::Mutation_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Mutation_Ptr_List_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Mutation_Ptr_List_Node_Traits", "get_next")] = lambda v: v['_mut_ptr_next']
boost.add_trivial_printer('MAC::Mutation_Ptr_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Read_Chunk_ITree_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_left")] = lambda v: v['_ce_l_child']
boost.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_right")] = lambda v: v['_ce_r_child']
boost.static_method[("MAC::detail::Read_Chunk_ITree_Node_Traits", "get_parent")] = lambda v: v['_ce_parent']
boost.add_trivial_printer('MAC::Read_Chunk_CE_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Read_Chunk_Set_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_left")] = lambda v: v['_re_l_child']
boost.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_right")] = lambda v: v['_re_r_child']
boost.static_method[("MAC::detail::Read_Chunk_Set_Node_Traits", "get_parent")] = lambda v: v['_re_parent']
boost.add_trivial_printer('MAC::Read_Chunk_RE_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Read_Chunk_Ptr_List_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Read_Chunk_Ptr_List_Node_Traits", "get_next")] = lambda v: v['_chunk_ptr_next']
boost.add_trivial_printer('MAC::Read_Chunk_Ptr_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Read_Entry_Set_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost.static_method[("MAC::detail::Read_Entry_Set_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost.add_trivial_printer('MAC::Read_Entry_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Contig_Entry_List_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Contig_Entry_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost.add_trivial_printer('MAC::Contig_Entry_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Hap_Entry_List_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Hap_Entry_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost.add_trivial_printer('MAC::Hap_Entry_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Hap_Hop_List_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Hap_Hop_List_Node_Traits", "get_next")] = lambda v: v['_next']
boost.add_trivial_printer('MAC::Hap_Hop_Cont', get_first_subclass)

boost.static_method[("MAC::detail::Hap_Hop_Set_Value_Traits", "to_value_ptr")] = identity
boost.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_left")] = lambda v: v['_l_child']
boost.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_right")] = lambda v: v['_r_child']
boost.static_method[("MAC::detail::Hap_Hop_Set_Node_Traits", "get_parent")] = lambda v: v['_parent']
boost.add_trivial_printer('MAC::Hap_Hop_Set', get_first_subclass)

mac_printer_gen = boost.Printer_Gen('MAC')
def add_mac_printer(Printer):
    mac_printer_gen.add(Printer)
    return Printer

@add_mac_printer
class Read_Chunk_Printer:
    printer_name = 'MAC::Read_Chunk'
    template_name = 'MAC::Read_Chunk'
    def __init__(self, v):
        self.v = v
        self.l = list()
        self.l.append(('read', '[%s,%s)' % (v['_r_start'], v['_r_start'] + v['_r_len'])))
        self.l.append(('contig', '[%s,%s)' % (v['_c_start'], v['_c_start'] + v['_c_len'])))
        self.l.append(('rc', str(boost.call_object_method(v, '_get_rc'))))
        self.l.append(('mut_ptr_cont', str(v['_mut_ptr_cont'])))
    def to_string(self):
        return None
    def children(self):
        return iter(self.l)

@add_mac_printer
class Read_Entry_Printer:
    printer_name = 'MAC::Read_Entry'
    template_name = 'MAC::Read_Entry'
    def __init__(self, v):
        self.v = v
        self.l = list()
        self.l.append(('name', v['_name']))
        self.l.append(('start', v['_start']))
        self.l.append(('len', v['_len']))
        self.l.append(('chunk_cont', str(v['_chunk_cont'])))
    def to_string(self):
        return None
    def children(self):
        return iter(self.l)

@add_mac_printer
class Mutation_Chunk_Adapter_Printer:
    printer_name = 'MAC::Mutation_Chunk_Adapter'
    template_name = 'MAC::Mutation_Chunk_Adapter'
    def __init__(self, v):
        self.v = v
        self.l = list()
        self.l.append(('mut_bptr', v['_mut_cbptr']))
        self.l.append(('chunk_bptr', v['_chunk_cbptr']))
    def to_string(self):
        return None
    def children(self):
        return iter(self.l)

gdb.printing.register_pretty_printer(None, mac_printer_gen)

end

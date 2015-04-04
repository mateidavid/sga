#!/usr/bin/env python
source itree.gdb
source factory.gdb
python
assert 'boost' in sys.modules

boost.static_method[('List_Node_Traits<Value>', 'get_next')] = lambda v: v['_list_next']
boost.static_method[('List_Value_Traits<Value>', 'to_value_ptr')] = lambda v: v

@boost.add_to_dict(boost.static_method, ('ITree_Node_Traits', 'get_parent'))
def f(ntt, node_rptr):
    return node_rptr['_parent']

@boost.add_to_dict(boost.static_method, ('ITree_Node_Traits', 'get_left'))
def f(ntt, node_rptr):
    return node_rptr['_l_child']

@boost.add_to_dict(boost.static_method, ('ITree_Node_Traits', 'get_right'))
def f(ntt, node_rptr):
    return node_rptr['_r_child']

@boost.add_to_dict(boost.static_method, ('ITree_Value_Traits', 'to_value_ptr'))
def f(vtt, node_rptr):
    return node_rptr

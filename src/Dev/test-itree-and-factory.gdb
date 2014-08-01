python
assert 'boost_print' in sys.modules, 'boost_print not found'

from boost_print import *
from boost_print.utils import _add_to_dict

@_add_to_dict(static_method, ('ITree_Node_Traits', 'get_parent'))
def f(ntt, node_rptr):
    return node_rptr['_parent']

@_add_to_dict(static_method, ('ITree_Node_Traits', 'get_left'))
def f(ntt, node_rptr):
    return node_rptr['_l_child']

@_add_to_dict(static_method, ('ITree_Node_Traits', 'get_right'))
def f(ntt, node_rptr):
    return node_rptr['_r_child']

@_add_to_dict(static_method, ('ITree_Value_Traits', 'to_value_ptr'))
def f(vtt, node_rptr):
    return node_rptr

#!/usr/bin/env python
python
assert 'boost' in sys.modules, 'boost not found'


# dict of max values
# key: type name
# value: max value (representing null)
#
max_val = dict()

def get_null_idn_val(t):
    idx_type = boost.get_basic_type(t.strip_typedefs().template_argument(1))
    idx_type_name = str(idx_type)
    if idx_type_name not in max_val:
        #p = boost.utils._new.invoke(boost.parse_and_eval('(%s *)0' % type_name))
        #res = p['_idx']
        #boost.utils._del.invoke(p)
        res = boost.parse_and_eval('(' + idx_type_name + ')-1')
        max_val[idx_type_name] = int(res)
    return max_val[idx_type_name]

# add null value checkers for bounded identifier and pointer
#
@boost.add_to_dict(boost.null_dict, 'bounded::detail::Identifier')
def __aux_factory_is_null_identifier(p):
    return int(p['_idx']) == get_null_idn_val(p.type)

@boost.add_to_dict(boost.null_dict, 'bounded::Pointer')
def __aux_factory_is_null_pointer(p):
    return is_null(p['_idn'])

# trivial identifier printer
#
def Identifier_Printer(v):
    if boost.is_null(v):
        return 'null'
    else:
        return str(v['_idx'])

def Pointer_Printer(p):
    # disabled, get_raw_ptr() CRASHES gdb when printing frame variables
    boost.message('printing: b_ptr(' + str(p['_idn']['_idx']) + ')')
    raw_ptr = boost.get_raw_ptr(p) ### CRASHES gdb
    boost.message('got raw_ptr = ' + str(raw_ptr))
    return str(p['_idn']) + ' (' + hex(int(raw_ptr)) + ')'

# trivial printer for bounded identifier and pointer objects
#
boost.add_trivial_printer('bounded::detail::Identifier', Identifier_Printer)
boost.add_trivial_printer('bounded::Pointer', lambda p: str(p['_idn']))

active_ptr_addr = dict()

def get_active_ptr(p):
    t = boost.get_basic_type(p.type)
    unqual_val_type = t.template_argument(0).unqualified()
    idx_type = t.template_argument(1)
    storage_type = boost.lookup_type('bounded::detail::Storage<' + str(unqual_val_type) + ',' + str(idx_type) + '>')
    if storage_type.tag in active_ptr_addr:
        return boost.parse_and_eval('*(' + storage_type.tag + '**)' + str(active_ptr_addr[storage_type.tag]))
    else:
        cmd_str = storage_type.tag + '::active_ptr()'
        try:
            return boost.parse_and_eval(cmd_str)
        except:
            boost.message('call failed: ' + cmd_str)
            boost.long_message(
                'get_active_ptr',
                '\n\tto bypass call:\n' +
                '\t- find static variable address with, e.g.: nm -C <executable> | grep ::_active_ptr\n' +
                '\t- add it to gdb with, e.g.: py active_ptr_addr["bounded::detail::Storage<MAC::Read_Entry, unsigned int>"]=0x0000000000951dd0')
            raise gdb.error

def factory_raw_ptr(p, idx=None):
    """
    Transform bounded pointer into raw pointer.
    If `idx` is given, it is taken to be the index to dereference, rather than the one in `p`.
    In this case, `p` still provides the bounded pointer type.
    """
    val_type = boost.get_basic_type(p.type.template_argument(0))
    null_idn_val = get_null_idn_val(p['_idn'].type)
    if idx == None:
        idx = boost.intptr(p['_idn']['_idx'])
    if idx == null_idn_val:
        return boost.parse_and_eval('(' + str(val_type) + ' *)0')
    # get active storage pointer
    active_storage_ptr = get_active_ptr(p)
    # complete reverse-engineer deque storage
    deque_ptr = active_storage_ptr['_cont'].address
    deque_impl_ptr = deque_ptr['_M_impl'].address
    deque_map = deque_impl_ptr['_M_map']
    elem_size = boost.get_basic_type(deque_ptr.type.target()).template_argument(0).sizeof
    page_size = max(int(512/elem_size), 1)
    node_start = deque_impl_ptr['_M_start']['_M_node']
    node_end = deque_impl_ptr['_M_finish']['_M_node']
    node_start_start = deque_impl_ptr['_M_start']['_M_first']
    node_start_crt = deque_impl_ptr['_M_start']['_M_cur']
    node_start_end = deque_impl_ptr['_M_start']['_M_last']
    node_end_start = deque_impl_ptr['_M_finish']['_M_first']
    node_end_crt = deque_impl_ptr['_M_finish']['_M_cur']
    node_end_end = deque_impl_ptr['_M_finish']['_M_last']
    # map is occupied as follows:
    #   start node is occupied from [node_start_crt, node_start_end)
    #   internal nodes are occupied with page_size elements
    #   end_node is occupied from [node_end_start, node_end_crt)
    n_elem = (boost.intptr(node_end - node_start - 1) * page_size
              + boost.intptr(node_start_end - node_start_crt)
              + boost.intptr(node_end_crt - node_end_start))
    if idx >= n_elem:
        return boost.parse_and_eval('(' + str(p.type.template_argument(0)) + ' *)0')
    if idx < boost.intptr(node_start_end - node_start_crt):
        vw_ptr = node_start_crt + idx
    else:
        idx2 = idx - boost.intptr(node_start_end - node_start_crt)
        k = int(idx2 / page_size)
        node_crt = node_start + 1 + k
        vw_ptr = node_crt.dereference() + (idx2 - k * page_size)
    return vw_ptr['_val'].address

    sz = boost.parse_and_eval(
        'bounded::detail::Storage< ' +
        str(boost.get_basic_type(p.type.template_argument(0))) + ', ' +
        str(p.type.template_argument(1)) + ' >'
        #+ '::active_ptr()->_cont.size()')
        + '::active_ptr()->ns_size()')
        #+ '::size()')
    if idx >= int(sz):
        return boost.parse_and_eval('(' + str(p.type.template_argument(0)) + ' *)0')
    #res = boost.parse_and_eval(
    #    '(' + str(p.type.template_argument(0)) + ' *)(& bounded::detail::Storage< ' +
    #    str(gdb.types.get_basic_type(p.type.template_argument(0))) + ', ' +
    #    str(p.type.template_argument(1)) + ' >::active_ptr()->_cont.at(' + str(idx) + '))')
    res = boost.parse_and_eval(
        '& bounded::detail::Storage< ' +
        str(boost.get_basic_type(p.type.template_argument(0))) + ', ' +
        str(p.type.template_argument(1)) + ' >::elem_at(' + str(idx) + ')')
    #print('got: ' + str(res), file=sys.stderr)
    return res

@boost.add_to_dict(boost.raw_ptr, 'bounded::Pointer')
def __aux_factory_raw_ptr(p):
    return factory_raw_ptr(p)

# convenience function for obtaining raw pointers from bounded identifiers and pointers
#
class raw_func(gdb.Function):
    def __init__(self):
        super(raw_func, self).__init__('raw')
    def invoke(self, p, idx=None):
        assert isinstance(p, gdb.Value)
        assert boost.template_name(p.type) == 'bounded::Pointer', '"v" not a bounded pointer'
        return factory_raw_ptr(p, idx)

raw = raw_func()

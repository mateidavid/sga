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

def factory_raw_ptr(p, idx=None):
    """
    Transform bounded pointer into raw pointer.
    If `idx` is given, it is taken to be the index to dereference, rather than the one in `p`.
    In this case, `p` still provides the bounded pointer type.
    """
    null_idn_val = get_null_idn_val(p['_idn'].type)
    if idx == None:
        idx = boost.intptr(p['_idn']['_idx'])
    if idx == null_idn_val:
        return boost.parse_and_eval('(' + str(p.type.template_argument(0)) + ' *)0')
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

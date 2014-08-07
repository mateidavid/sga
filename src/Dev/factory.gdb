python
assert 'boost_print' in sys.modules, 'boost_print not found'

#from boost_print import *
#from boost_print.utils import _add_to_dict, _new, _del

# dict of null identifiers
# key: type name
# value: _idx value representing null
#
null_idn_val = dict()

def get_null_idn_val(t):
    type_name = str(t.strip_typedefs())
    if type_name not in null_idn_val:
        #gdb.execute('set $_tmp = (%s *)malloc(sizeof(%s))' % (type_name, type_name))
        #gdb.execute('call $_tmp->Identifier()')
        #res = gdb.parse_and_eval('$_tmp._idx')
        #gdb.execute('call free($_tmp)')
        p = boost_print.utils._new.invoke(boost_print.parse_and_eval('(%s *)0' % type_name))
        res = p['_idx']
        boost_print.utils._del.invoke(p)
        null_idn_val[type_name] = int(res)
    return null_idn_val[type_name]

# add null value checkers for bounded identifier and pointer
#
@boost_print.utils._add_to_dict(boost_print.null_dict, 'bounded::detail::Identifier')
def __aux_factory_is_null_identifier(p):
    return int(p['_idx']) == get_null_idn_val(p.type)

@boost_print.utils._add_to_dict(boost_print.null_dict, 'bounded::detail::Pointer')
def __aux_factory_is_null_pointer(p):
    return boost_print.is_null(p['_idn'])

# trivial identifier printer
#
def Identifier_Printer(v):
    if boost_print.is_null(v):
        return 'null'
    else:
        return str(v['_idx'])

def Pointer_Printer(p):
    # disabled, get_raw_ptr() CRASHES gdb when printing frame variables
    print('printing: b_ptr(' + str(p['_idn']['_idx']) + ')', file=sys.stderr)
    raw_ptr = boost_print.get_raw_ptr(p) ### CRASHES gdb
    print('got raw_ptr = ' + str(raw_ptr), file=sys.stderr)
    return str(p['_idn']) + ' (' + hex(int(raw_ptr)) + ')'

# trivial printer for bounded identifier and pointer objects
#
boost_print.add_trivial_printer('bounded::detail::Identifier<.*>', Identifier_Printer)
boost_print.add_trivial_printer('bounded::detail::Pointer<.*>', lambda p: str(p['_idn']))

@boost_print.utils._add_to_dict(boost_print.raw_ptr, 'bounded::detail::Pointer')
def __aux_factory_raw_ptr(p):
    if boost_print.is_null(p):
        return boost_print.parse_and_eval(
            '(' + str(p.type.template_argument(0).strip_typedefs()) + ' *)0')

    idn = p['_idn']
    idx = idn['_idx']
    #print('raw_ptr(' + str(idx) + ')', file=sys.stderr)
    idn_str = boost_print.to_eval(idn, '$_raw_bptr_idn')
    #cmd = ('& bounded::detail::Storage<'
    #       + str(p.type.template_argument(0)) + ', ' + str(p.type.template_argument(1))
    #       + '>::elem_at(' + idn_str + ')')
    sz = boost_print.parse_and_eval(
        'bounded::detail::Storage< ' +
        str(p.type.template_argument(0)) + ', ' + str(p.type.template_argument(1)) +
        ' >::_active_ptr->_cont.size()')
    if idx >= int(sz):
        return boost_print.parse_and_eval(
            '(' + str(p.type.template_argument(0)) + ' *)0')
    res = boost_print.parse_and_eval(
        '(' + str(p.type.template_argument(0)) + ' *)(& bounded::detail::Storage< ' +
        str(p.type.template_argument(0)) + ', ' + str(p.type.template_argument(1)) +
        ' >::_active_ptr->_cont.at(' + str(idx) + '))')
    #print('got: ' + str(res), file=sys.stderr)
    return res


# convenience function for obtaining raw pointers from bounded identifiers and pointers
#
class raw_func(gdb.Function):
    def __init__(self):
        super(raw_func, self).__init__('raw')
    def invoke(self, p):
        assert isinstance(p, gdb.Value), '"v" not a gdb.Value'
        assert boost_print.template_name(p.type) == 'bounded::detail::Pointer', '"v" not a bounded pointer'
        return boost_print.get_raw_ptr(p)

raw = raw_func()

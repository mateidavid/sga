python
assert 'boost_print' in sys.modules, 'boost_print not found'

from boost_print import *
from boost_print.utils import _add_to_dict, _new, _del

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
        p = _new.invoke(gdb.parse_and_eval('(%s *)0' % type_name))
        res = p['_idx']
        _del.invoke(p)
        null_idn_val[type_name] = int(res)
    return null_idn_val[type_name]

# add null value checkers for bounded identifier and pointer
#
@_add_to_dict(null_dict, 'bounded::detail::Identifier')
def f(p):
    return int(p['_idx']) == get_null_idn_val(p.type)

@_add_to_dict(null_dict, 'bounded::detail::Pointer')
def f(p):
    return is_null(p['_idn'])

# trivial identifier printer
#
def Identifier_Printer(v):
    if is_null(v):
        return 'null'
    else:
        return str(v['_idx'])

# trivial printer for bounded identifier and pointer objects
#
add_trivial_printer('bounded::detail::Identifier<.*>', Identifier_Printer)
add_trivial_printer('bounded::detail::Pointer<.*>',
                    lambda p: str(p['_idn']) + ' (' + hex(int(get_raw_ptr(p))) + ')')

@_add_to_dict(raw_ptr, 'bounded::detail::Pointer')
def f(p):
    if is_null(p):
        return gdb.parse_and_eval(
            '(' + str(p.type.template_argument(0).strip_typedefs()) + ' *)0')

    idn = p['_idn']
    idn_str = to_eval(idn, '$_idn')
    return gdb.parse_and_eval(
        '& bounded::detail::Storage<'
        + str(p.type.template_argument(0)) + ', ' + str(p.type.template_argument(1))
        + '>::elem_at(' + idn_str + ')')


# convenience function for obtaining raw pointers from bounded identifiers and pointers
#
class raw_func(gdb.Function):
    def __init__(self):
        super(raw_func, self).__init__('raw')
    def invoke(self, p):
        assert isinstance(p, gdb.Value), '"v" not a gdb.Value'
        assert template_name(p.type) == 'bounded::detail::Pointer', '"v" not a bounded pointer'
        return get_raw_ptr(p)

raw = raw_func()

#ifndef __FOR_EACH_ADVANCE_HPP
#define __FOR_EACH_ADVANCE_HPP


/**
 * Apply function to elements in range, but advance iterator first.
 * With this, it is safe to apply a function which might indirectly remove
 * the pointed element from a list.
 * @param first Start of range.
 * @param last Enf of range.
 * @param f Unary function taking an element as argument.
 * @return f
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function for_each_advance(Input_Iterator first, Input_Iterator last, Unary_Function fn)
{
    while (first != last)
    {
        fn(*(first++));
    }
    return fn;
}

/**
 * Apply function to iterators in range.
 * @param first Start of range.
 * @param last Enf of range.
 * @param f Unary function taking an iterator as argument.
 * @return f
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function for_each_it(Input_Iterator first, Input_Iterator last, Unary_Function fn)
{
    while (first != last)
    {
        fn(first);
        ++first;
    }
    return fn;
}

/**
 * Apply function to iterators in range, but advance iterator first.
 * With this, it is safe to remove elements from a list on the go.
 * @param first Start of range.
 * @param last Enf of range.
 * @param f Unary function taking an element as argument.
 * @return f
 */
template < class Input_Iterator, class Unary_Function >
Unary_Function for_each_it_advance(Input_Iterator first, Input_Iterator last, Unary_Function fn)
{
    while (first != last)
    {
        fn(first++);
    }
    return fn;
}


#endif

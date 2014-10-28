//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __BITS_HPP
#define __BITS_HPP


namespace bitmask
{

/* Clear the bits set in the bitmask. */
template < typename T >
void reset(T& dest, T mask)
{
    dest &= ~mask;
}

/* Set the bits set in the bitmask. */
template < typename T >
void set(T& dest, T mask)
{
    dest |= mask;
}

/* Set or reset the bits set in the bitmask, depending on val. */
template < typename T >
void assign(T& dest, T mask, bool val)
{
    if (val)
    {
        set(dest, mask);
    }
    else
    {
        reset(dest, mask);
    }
}

/* Check if any of the bits set in the bitmask is set */
template < typename T >
bool any(const T& dest, T mask)
{
    return dest & mask;
}

/* Check if all of the bits set in the bitmask are set */
template < typename T >
bool all(const T& dest, T mask)
{
    return (dest & mask) == mask;
}

}


#endif

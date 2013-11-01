//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
// Inspired by:
// http://stackoverflow.com/questions/13842676/c-custom-output-stream-with-indentation
//-----------------------------------------------

#ifndef __INDENT_HPP
#define __INDENT_HPP

#include <iostream>


namespace indent
{
    /** Manipulator that increases indent level. */
    std::ios_base& inc(std::ios_base&);

    /** Manipulator that decreases indent level. */
    std::ios_base& dec(std::ios_base&);

    /** Manipulator that inserts current indent level. */
    std::ostream& tab(std::ostream& os);

    /** Manipulator that inserts new line and indents. */
    std::ostream& nl(std::ostream& os);

    /** Manipulator that inserts std::endl (which flushes) and indents. */
    std::ostream& endl(std::ostream& os);
}


#endif

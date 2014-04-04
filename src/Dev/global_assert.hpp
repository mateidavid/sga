//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#ifndef __GLOBAL_ASSERT_HPP
#define __GLOBAL_ASSERT_HPP

#include <iostream>


#define BOOST_ENABLE_ASSERT_HANDLER
#define BOOST_ASSERT_MSG_OSTREAM std::cerr
#include <boost/assert.hpp>
#define ASSERT BOOST_ASSERT
#define ASSERT_MSG BOOST_ASSERT_MSG


#endif

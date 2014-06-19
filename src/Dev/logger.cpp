//-----------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Matei David (mdavid@oicr.on.ca)
// Released under the GPL license
//-----------------------------------------------

#include "logger.hpp"


namespace logger
{

std::map< std::string, level > log::_facility_level;
level log::_default_level = error;

} // namespace logger

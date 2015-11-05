#ifndef structmap_hh
#define structmap_hh

#ifndef structmap_pp_test
#include <iostream>
#include <string>
#include <stdexcept>
#endif

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#define STRUCTMAP_STR_CASE(r, enum_name, elem) \
  if (!str.compare( BOOST_PP_STRINGIZE(elem) )) return elem; else

#define structmap(type, map_name, values) \
  struct map_name { \
    type BOOST_PP_SEQ_ENUM(values); \
    type& operator[](const std::string& str) { \
      BOOST_PP_SEQ_FOR_EACH( STRUCTMAP_STR_CASE, nil, values ) \
      throw std::runtime_error( \
        "structmap " BOOST_PP_STRINGIZE(map_name) \
        " does not map \""+str+"\"" \
      ); \
    } \
  }

#endif

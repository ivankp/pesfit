#ifndef senum_hh
#define senum_hh

#ifndef senum_pp_test
#include <iostream>
#include <string>
#include <stdexcept>
#endif

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/comparison/greater.hpp>

#define SENUM_EN_CASE(r, data, elem) \
  case elem: return BOOST_PP_STRINGIZE(elem);

#define SENUM_STR_CASE(r, enum_name, elem) \
  if (!str.compare( BOOST_PP_STRINGIZE(elem) )) return elem; else

#define SENUM_CS(r, enum_name, elem) BOOST_PP_IF(BOOST_PP_GREATER(r,2),", ",) BOOST_PP_STRINGIZE(elem)

#define senum(enum_name, values) \
  struct enum_name { \
    enum type { BOOST_PP_SEQ_ENUM(values) }; \
    enum { _nelem = BOOST_PP_SEQ_SIZE(values) }; \
    static type _enum(const std::string& str) { \
      BOOST_PP_SEQ_FOR_EACH( SENUM_STR_CASE, nil, values ) \
      throw std::runtime_error( \
        "invalid value \""+str+"\" specified for stringized enum " \
        BOOST_PP_STRINGIZE(enum_name) \
      ); \
    } \
    static const char* _str(type en) { \
      switch (en) { \
        BOOST_PP_SEQ_FOR_EACH( SENUM_EN_CASE, nil, values ) \
      } \
      return 0; \
    } \
    static const char* _str_all() { \
      return BOOST_PP_SEQ_FOR_EACH( SENUM_CS, nil, values ); \
    } \
  }; \
\
std::istream& operator>> (std::istream& is, enum_name::type& x) { \
  std::string str; \
  is >> str; \
  x = enum_name::_enum(str); \
  return is; \
} \
\
std::ostream& operator<< (std::ostream& os, const enum_name::type& x) { \
  os << enum_name::_str(x); \
  return os; \
}

#endif

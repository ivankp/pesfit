#ifndef regex_hh
#define regex_hh

#define GCC_VERSION (__GNUC__ * 10000 \
                              + __GNUC_MINOR__ * 100 \
                              + __GNUC_PATCHLEVEL__)

#if GCC_VERSION < 40900 // Test for GCC < 4.9.0
#include <boost/regex.hpp>
using boost::regex;
using boost::smatch;
using boost::regex_match;
#define regex_icase boost::regex::icase
#else
#include <regex>
using std::regex;
using std::smatch;
using std::regex_match;
using std::regex_constants::icase;
#define regex_icase std::regex_constants::icase
#endif

#endif

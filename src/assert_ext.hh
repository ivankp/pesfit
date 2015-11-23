#ifndef assert_ext_hh
#define assert_ext_hh

#include <vector>
#include <initializer_list>
#include <sstream>
#include <stdexcept>

bool check_ext(
  const std::string& fname,
  const std::initializer_list<const char*>& exts
) {
  const std::string fext = fname.substr(fname.rfind('.')+1);
  for (const auto& ext : exts)
    if (!fext.compare(ext)) return true;
  return false;
}

bool check_ext(
  const std::vector<std::string>& fnames,
  const std::initializer_list<const char*>& exts
) {
  for (const auto& fname : fnames)
    if (!check_ext(fname, exts)) return false;
  return true;
}

void assert_ext(
  const std::string& fname,
  const std::initializer_list<const char*>& exts
) {
  if (!check_ext(fname, exts)) {
    std::stringstream ss;
    ss << "Unexpected extension in filename: "<<fname<<'\n';
    ss << "Allowed extensions are:";
    for (const auto& ext : exts) ss << ' ' << ext;
    throw std::runtime_error(ss.str());
  }
}

void assert_ext(
  const std::vector<std::string>& fnames,
  const std::initializer_list<const char*>& exts
) {
  for (const auto& fname : fnames)
    assert_ext(fname, exts);
}

void assert_ext(
  const char* const * fnames, int n,
  const std::initializer_list<const char*>& exts
) {
  for (int i=0; i<n; ++i)
    assert_ext(fnames[i], exts);
}

#endif

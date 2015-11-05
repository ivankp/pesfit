#ifndef workspace_hh
#define workspace_hh

#include <string>
#include <exception>
class root_get_error: public std::exception {
  std::string what_arg;
public:
  explicit root_get_error(const std::string& what_arg): what_arg(what_arg) { }
  explicit root_get_error(const char* what_arg): what_arg(what_arg) { }
  virtual const char* what() const noexcept { return what_arg.c_str(); }
};

class TObject;
class TDirectory;

template<typename T>
inline T* get(TDirectory* d, const char* name) {
  TObject *obj = d->Get(name);
  if (!obj) {
    throw root_get_error( std::string("No object ")
      +name+" in "+d->GetName()
    );
  } else if (obj->InheritsFrom(T::Class())) {
    return static_cast<T*>(obj);
  } else {
    throw root_get_error( std::string(obj->ClassName())
      +' '+obj->GetName()+" does not inherit from "
      +T::Class()->GetName()
    );
  }
}

#endif

#ifndef root_canvas_table_hh
#define root_canvas_table_hh

#include <vector>
#include <sstream>

class TText;
class TPaveText;

class table {
  float x1, dy1, x2, dy2;
  std::vector<std::vector<TPaveText*>> cells;
public:
  table(int nrows, int ncols, float x1=0., float y1=0., float x2=1., float y2=1.)
  ~table();

  TPaveText* Get(size_t i, size_t j);
  TText* operator()(size_t i, size_t j, const char* text);

  struct stream {
    std::stringstream ss;
    TPaveText* pave;
    stream(TPaveText* pave): pave(pave) { }
    ~stream()
    template<typename T>
    inline std::istream& operator<<(T&& x) { return ss << x; }
  };
  stream operator()(size_t i, size_t j);

  void Draw(bool outer=false, bool inner=true);
};

#endif

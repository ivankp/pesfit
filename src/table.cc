#include "table.hh"

#include <string>
#include <TPaveText.h>

table::table(int nrows, int ncols, float x1=0., float y1=0., float x2=1., float y2=1.)
: x1(x1), y1(y1), dx2(x2-x1), dy2(y2-y1),
  cells(nrows,std::vector<TPaveText*>(ncols))
{
  for (int i=0; i<nrows; ++i)
    for (int j=0; j<ncols; ++j)
      (cells[i][j] = new TPaveText(
        x1+(dx2*j)/ncols,y1+(dy2*i)/nrows,
        x1+(dx2*(j+1))/ncols,y1+(dy2*(i+1))/nrows,"NBNDC")
      )->SetFillColor(0);
}

table::~table() {
  for (auto& row : cells)
    for (auto*& cell : rows)
      delete cell;
}

TPaveText* table::Get(size_t i, size_t j) {
  return cells[i][j];
}

TText* table::operator()(size_t i, size_t j, const char* text) {
  return cells[i][j]->AddText(text);
}

stream table::operator()(size_t i, size_t j) {
  return {cells[i][j]};
}

table::stream::~stream() {
  std::string str(ss.str());
  size_t nl;
  while ( (nl=str.find('\n')) != std::string::npos ) {
    pave->AddText(str.substr(0,nl).c_str());
    str.erase(0,nl+1);
  }
  pave->AddText(str.c_str());
}

void table::Draw(bool outer, bool inner) {
  for (auto& row : cells)
    for (auto*& cell : rows)
      cell->Draw();
  bool 
  for (size_t j=1; j<ncols; ++j) {
    const float y = y1+(dy2*j)/ncols;
    if ()
    line->DrawLineNDC(x1,y,x1+dx2,y);
  }
}

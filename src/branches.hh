// Developed by Ivan Pogrebnyak, MSU

#ifndef ttree_branches_hh
#define ttree_branches_hh

#include <TTree.h>

template<typename T>
inline void branches_impl(TTree* tree, const char* name, T* add) {
  tree->SetBranchAddress(name, add);
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus(name,1);
}

template<typename T, typename... TT>
inline void branches_impl(
  TTree* tree, const char* name, T* add, TT... bb
) {
  tree->SetBranchAddress(name, add);
  branches(tree,bb...);
  tree->SetBranchStatus(name,1);
}

template<typename... TT>
inline void branches(TTree* tree, TT... bb) {
  tree->SetBranchStatus("*",1);
  branches_impl(tree,bb...);
}

#endif

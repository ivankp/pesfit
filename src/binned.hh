#ifndef binned_hh
#define binned_hh

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

template<typename T, typename Bin=double, typename Compare=std::less<Bin>>
class binned {
public:
  typedef Bin     bin_type;
  typedef T       value_type;
  typedef Compare bin_compare;

public:
  std::vector<bin_type>   bins;
  std::vector<value_type> vals;

public:
  binned(const std::vector<bin_type>& bins)
  : bins(bins), vals(bins.size()+1)
  {
    if (!is_sorted(bins.begin(),bins.end()))
      throw std::runtime_error("Unsorted bins vector passed to binned container");
  }

  binned(const binned<value_type,bin_type,bin_compare>& other, bool blank)
  : bins(other.bins), vals(bins.size()+1)
  {
    if (!blank)
      vals.assign(other.vals.begin(),other.vals.end());
  }

  template<typename Min, typename Step>
  binned(size_t nbins, const Min& min, const Step& step)
  : bins(nbins), vals(nbins+1)
  {
    bins[0] = min;
    for (size_t i=1; i<nbins; ++i)
      bins[i] = bins[i-1] + step;
  }

  value_type& operator[](const bin_type& x) {
    if ( bin_compare()(x,bins.front()) ) {
      return vals.front();
    } else if ( !bin_compare()(x,bins.back()) ) {
      return vals.back();
    } else {
      auto bin = upper_bound(bins.begin(),bins.end(),x,bin_compare());
      return *(vals.begin()+((bin-bins.begin())));
    }
  }

  auto begin(bool underflow=false) -> typename decltype(vals)::iterator {
    if (underflow) return vals.begin();
    else return vals.begin()+1;
  }

  auto end(bool overflow=false) -> typename decltype(vals)::iterator {
    return vals.begin()+(overflow?bins.size()+1:bins.size());
  }

  inline size_t nbins() const noexcept {
    return bins.size()-1;
  }

  inline value_type& at(size_t i) {
    return vals.at(i);
  }

  size_t bin_index(const bin_type& x) const {
    if ( bin_compare()(x,bins.front()) ) {
      return 0;
    } else if ( !bin_compare()(x,bins.back()) ) {
      return bins.size();
    } else {
      auto bin = upper_bound(bins.begin(),bins.end(),x,bin_compare());
      return bin-bins.begin();
    }
  }

  const bin_type& left_edge(size_t i) const {
    if (i==0)
      throw std::runtime_error("Bin 0 is underflow");
    else if (i > bins.size())
      throw std::runtime_error("Bin "+std::to_string(bins.size())+" is overflow");
    return bins[i-1];
  }
  const bin_type& right_edge(size_t i) const {
    if (i >= bins.size())
      throw std::runtime_error("Bin "+std::to_string(bins.size())+" is overflow");
    return bins[i];
  }

};

#endif

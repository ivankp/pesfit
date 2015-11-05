#ifndef seqmap_hh
#define seqmap_hh

template<typename T, typename Key=std::string>
class seqmap: public std::vector<std::pair<Key,T>> {
public:
  T& operator[](const Key& key) {
    auto it = std::find_if(this->begin(),this->end(),
      [&key](const std::pair<Key,T>& p){ return p.first==key; });
    if (it==this->end()) {
      this->emplace_back(key,T());
      return this->back().second;
    } else return it->second;
  }
};

#endif

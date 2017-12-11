#ifndef RANDOMACCESSSEQUENCE_H
#define RANDOMACCESSSEQUENCE_H

#include <functional>

/** 
    A generic random access sequence type whose values are stored externally
 */
template < typename _RandomAccessInputIterator >
class RandomAccessSequence {
  _RandomAccessInputIterator _begin;
  _RandomAccessInputIterator _end;
  size_t _len;
public:
  typedef typename std::iterator_traits<_RandomAccessInputIterator>::value_type ElemTy;
  RandomAccessSequence() {}
  RandomAccessSequence(_RandomAccessInputIterator b,
                       _RandomAccessInputIterator e):
    _begin(b), _end(e), _len(e-b) {}
  ~RandomAccessSequence() {}
  inline _RandomAccessInputIterator begin() { return _begin; }
  inline _RandomAccessInputIterator end() { return _end; }
  inline ElemTy pop_front() { ElemTy tmp = *_begin; ++_begin; --_len; return tmp; }
  inline ElemTy pop_back() {  ElemTy tmp = *(_end-1); --_end; --_len;  return tmp; }
  inline size_t size() const { return _len; }
  inline ElemTy operator[] (size_t index) const {assert(index < _len); return *(_begin + index); }

  template <typename Equivalent = std::equal_to<void>>
      unsigned find(ElemTy elem, Equivalent cmp = Equivalent()) {
    unsigned index = 0;
    for (_RandomAccessInputIterator i = _begin; i < _end; ++i, ++index)
      if (cmp(*i, elem))
        return index;
    return -1;
  }
  void split(size_t index, RandomAccessSequence &left, RandomAccessSequence &right) {
    left = RandomAccessSequence(_begin, _begin + index);
    right = RandomAccessSequence(_begin + index, _end);
  }
};
 
#endif

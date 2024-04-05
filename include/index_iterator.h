#ifndef INDEX_ITERATOR
#define INDEX_ITERATOR


#include "index_type.h"





class IndexIterator {
public:
   IndexIterator() :
      begin(0), loc(0), size(0) {};
   IndexIterator(index_type* b, index_type s) :
      begin(b), loc(b), size(s) {};
   IndexIterator(const IndexIterator& ii) :
      begin(ii.begin), loc(ii.loc), size(ii.size) {};

   index_type* begin;
   index_type* loc;
   index_type size;

   inline bool isValid() {
      return (loc-begin < size);
   }

   //prefix
   inline IndexIterator& operator++();
   //postfix
   inline IndexIterator operator++(int);
};


inline IndexIterator& IndexIterator::operator++() {
   if (loc-begin < size) {
      loc++;
   }
   return *this;
}

inline IndexIterator IndexIterator::operator++(int) {
   IndexIterator ret = *this;
   ++(*this);
   return ret;
}



#endif

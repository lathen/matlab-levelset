/*************************************************************************************************
 *
 * Copyright (c) 2008
 * Gunnar Läthén (gunnar.lathen@itn.liu.se)
 * Linköping University, Sweden.
 *
 *************************************************************************************************
 * Contributors:  
 *                  1) Gunnar Läthén (gunnar.lathen@itn.liu.se)
 *************************************************************************************************
 *
 * This class implements a simple heap data structure.
 *
 *************************************************************************************************/

#ifndef _HEAP_H
#define _HEAP_H

#include <vector>
#include <iostream>
#include <limits>

template <typename DataType>
class Heap
{
public :

  Heap();
  ~Heap();

  class Heapable
  {
  public :
    Heapable() : cost(0), position((std::numeric_limits<unsigned int>::max)()) { }
    Heapable(DataType val) : cost(val), position((std::numeric_limits<unsigned int>::max)()) { }

    DataType cost;
    unsigned int position;

    bool operator < (const Heapable & h) const {
      return this->cost < h.cost;
    }
  };

  void push(Heapable * h);

  Heapable * peek();

  Heapable * pop();

  Heapable * remove(Heapable * h);
  
  void reserve(unsigned int n);

  inline unsigned int size() { return mNodes.size()-1; }

  inline bool isEmpty() { return size() == 0; }

  void update(Heapable * h);

  void print(std::ostream & os);

protected :

  inline unsigned int parent(unsigned int i) { return i/2; }
  inline unsigned int leftChild(unsigned int i) { return 2*i; }
  inline unsigned int rightChild(unsigned int i) { return 2*i+1; }

  void percolateUp(unsigned int hole);
  void percolateDown(unsigned int hole);

  std::vector<Heapable *> mNodes;
};

#include "Heap_Impl.h"

#endif

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
 * This is the implementation file for the heap data structure.
 *
 *************************************************************************************************/

#include <cassert>

template <typename DataType>
Heap<DataType>::Heap()
{
  // Add a sentinel (dummy node) for convinience
  mNodes.push_back(new Heapable(-(std::numeric_limits<DataType>::max)()));
}

template <typename DataType>
Heap<DataType>::~Heap()
{
  // delete sentinel node
  delete mNodes.front();
}


template <typename DataType>
void Heap<DataType>::push(Heapable * h)
{
  mNodes.push_back(h);
  percolateUp(mNodes.size()-1);
}


template <typename DataType>
typename Heap<DataType>::Heapable * Heap<DataType>::peek()
{
  return (mNodes.size() > 1 ? mNodes[1] : NULL);
}


template <typename DataType>
typename Heap<DataType>::Heapable * Heap<DataType>::pop()
{
  if (mNodes.size() == 1)  return NULL;
  return remove(mNodes[1]);
}


template <typename DataType>
typename Heap<DataType>::Heapable * Heap<DataType>::remove(Heapable * h)
{
  assert(h->position != (std::numeric_limits<unsigned int>::max)() && h->position < mNodes.size());

  unsigned int hole = h->position;
  mNodes[hole] = mNodes[mNodes.size()-1];
  mNodes.pop_back();
  h->position = (std::numeric_limits<unsigned int>::max)();

  if (hole == mNodes.size()) return h;

  if (*mNodes[hole] < *h)
    percolateUp(hole);
  else
    percolateDown(hole);

  return h;
}


template <typename DataType>
void Heap<DataType>::reserve(unsigned int n)
{
  mNodes.reserve(n);
}


template <typename DataType>
void Heap<DataType>::update(Heapable * h)
{
  assert(h->position != (std::numeric_limits<unsigned int>::max)() && h->position < mNodes.size());

  if (*h < *mNodes[parent(h->position)])
    percolateUp(h->position);
  else
    percolateDown(h->position);
}


template <typename DataType>
void Heap<DataType>::print(std::ostream & os)
{
  typename std::vector<Heapable *>::const_iterator iter = mNodes.begin();
  typename std::vector<Heapable *>::const_iterator iend = mNodes.end();
  ++iter; // skip the sentinel
  while (iter != iend) {
    os << (*iter)->cost << "(" << (*iter)->position << ") ";
    ++iter;
  }
  os << std::endl;
}


template <typename DataType>
void Heap<DataType>::percolateUp(unsigned int hole)
{
  Heapable * start = mNodes[hole];
  while (*start < *mNodes[ parent(hole) ]) {
    mNodes[hole] = mNodes[ parent(hole) ];
    mNodes[hole]->position = hole;
    hole = parent(hole);
  }
  mNodes[hole] = start;
  start->position = hole;
}

template <typename DataType>
void Heap<DataType>::percolateDown(unsigned int hole)
{
  Heapable * start = mNodes[hole];
  unsigned int currentSize = mNodes.size();

  while (leftChild(hole) < currentSize) {
    unsigned int left = leftChild(hole);
    unsigned int right = rightChild(hole);
    unsigned int child = left;
    if (right < currentSize && *mNodes[right] < *mNodes[left])
      child = right;

    if (*mNodes[child] < *start) {
      mNodes[hole] = mNodes[child];
      mNodes[hole]->position = hole;
    }
    else break;

    hole = child;
  }
  mNodes[hole] = start;
  start->position = hole;
}

/*
 * small_object_storage.hpp
 *
 *  Created on: 18 Nov 2020
 *  Author: Holger Schmitz (holger@notjustphysics.com)
 */

#ifndef HUERTO_STORAGE_SMALL_OBJECT_STORAGE_HPP_
#define HUERTO_STORAGE_SMALL_OBJECT_STORAGE_HPP_

#include <list>
#include <memory>
#include <fstream>

struct DefaultSmallObjectStorageTraits {
    static const size_t chunkSize = 1000000;
};

/**
 * @brief A dynamic storage for small(ish) objects of equal size.
 * 
 * The objects are stored in chunks of a specified size.
 * The chunks are stored in a STL list.
 * 
 * Insertion of elements is performed by writing the element into the next space of
 * a chunk that contains free spaces. Internally a reference of such a chunk is held.
 * This means, most insertions are performed in constant time.
 * 
 * Provides an STL style iterator and `begin()` and `end()` methods to iterate over all 
 * elements.
 */
template<typename T, class Traits = DefaultSmallObjectStorageTraits>
class SmallObjectStorage {
  public:
    typedef T value_type;
  private:
    struct DataChunk {
        T *data;
        size_t count;
        DataChunk() = delete;
        DataChunk(T *data);
        DataChunk(const DataChunk &chunk);
        T &addItem();
    };

    typedef std::list<DataChunk> ChunkList;
    typedef typename ChunkList::iterator ChunkIterator;
    typedef typename ChunkList::const_iterator ChunkConstIterator;
    ChunkList chunks;
    ChunkIterator freeChunk;
  public:

    SmallObjectStorage() {
//      chunks.push_front(DataChunk(new T[Traits::chunkSize]));
      freeChunk = chunks.end();
    }
    
    SmallObjectStorage(const SmallObjectStorage&) = delete;

    ~SmallObjectStorage();

    class iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
      private:
        friend class SmallObjectStorage<T, Traits>;
        ChunkIterator chunkIter;
        size_t pos;

        iterator(ChunkIterator chunkIter_, long pos_ = 0) :
          chunkIter(chunkIter_), pos(pos_) {}
      public:
        iterator() : pos(0) {}

        iterator(const iterator &it) : chunkIter(it.chunkIter), pos(it.pos) {}

        iterator& operator++() {
          if (++pos >= chunkIter->count) {
            ++chunkIter;
            pos = 0;
          }
          return *this;
        }

        iterator operator++(int) {
          iterator tmp(*this);
          operator++();
          return tmp;
        }

        iterator& operator--() {
          if (--pos < 0) {
            --chunkIter;
            pos = chunkIter->count - 1;
          }
          return *this;
        }

        iterator operator--(int) {
          iterator tmp(*this);
          operator--();
          return tmp;
        }

        bool operator==(const iterator& rhs) const {
          return (chunkIter == rhs.chunkIter) && (pos == rhs.pos);
        }

        bool operator!=(const iterator& rhs) const {
          return (chunkIter != rhs.chunkIter) || (pos != rhs.pos);
        }

        T& operator*() {
          return chunkIter->data[pos];
        }

        T* operator->() {
          return &(chunkIter->data[pos]);
        }
    };

    iterator begin() {
      ChunkIterator b = chunks.begin();
      while (b->count == 0 && b != chunks.end()) {
        ++b;
      }
      return iterator(b);
    }

    iterator end() {
      return iterator(chunks.end());
    }

    T &addItem();

    /**
     * Remove an item from the storage.
     *
     * After deletion the iterator will point to the position after the deleted item.
     * This is in line with STL behaviour.
     *
     * Iterators pointing to items before the deleted item remain valid.
     * Iterators pointing to items after the deleted item become invalid.
     */
    iterator removeItem(const iterator&);

    long getCount() const;

    /**
     * Empty the storage
     */
    void clear();
};

//=================================================================
//=============== SmallObjectStorage ==============================
//=================================================================

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::DataChunk::DataChunk(T *data) : data(data), count(0) {
}

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::DataChunk::DataChunk(const DataChunk &chunk) : data(chunk.data), count(chunk.count) {
}

template<typename T, class Traits>
T &SmallObjectStorage<T, Traits>::DataChunk::addItem() {
  return data[count++];
}

template<typename T, class Traits>
T &SmallObjectStorage<T, Traits>::addItem() {
  // find a free chunk or create one
  if ((freeChunk == chunks.end()) || (Traits::chunkSize == freeChunk->count)) {
    freeChunk = chunks.begin();

    while((freeChunk != chunks.end()) && (Traits::chunkSize <= freeChunk->count)) {
      ++freeChunk;
    }

    if (freeChunk == chunks.end()) {
      freeChunk = chunks.insert(chunks.begin(), DataChunk(new T[Traits::chunkSize]));
    }
  }
  return freeChunk->addItem();
}

template<typename T, class Traits>
typename SmallObjectStorage<T, Traits>::iterator SmallObjectStorage<T, Traits>::removeItem(const iterator &it_) {
  iterator it = it_;
  DataChunk &dchunk = *(it.chunkIter);
  T *data = dchunk.data;
  long pos = it.pos;
  long last = dchunk.count - 1;

  if (last>pos) std::swap(data[pos], data[last]);

  --(dchunk.count);
  // if no more items in the chunk then delete the chunk
  if (0 == dchunk.count && chunks.size() > 1) {
    std::cerr << "SmallObjectStorage deleting data" << dchunk.data << std::endl;
    delete[] dchunk.data;
    it.chunkIter = chunks.erase(it.chunkIter);
    it.pos = 0;
    freeChunk = chunks.begin();
  } else if (last==pos) {
    ++it.chunkIter;
    it.pos = 0;
  }

  // for performance, see if there are more free spaces in this chunk than in the
  // freeChunk, and update freeChunk accordingly
  //  if ((it.chunkIter != chunks.end()) && (it.chunkIter->count < freeChunk->count))
  //  {
  //    freeChunk = it.chunkIter;
  //  }
  return it;
}

template<typename T, class Traits>
long SmallObjectStorage<T, Traits>::getCount() const
{
  long count = 0;
  for (ChunkConstIterator b = chunks.begin(); b!=chunks.end(); ++b)
    count += b->count;
  return count;
}

template<typename T, class Traits>
void SmallObjectStorage<T, Traits>::clear() {
    for (ChunkIterator b = chunks.begin(); b!=chunks.end(); ++b) {
      delete[] b->data;
    }
    chunks.clear();
    freeChunk = chunks.end();
}

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::~SmallObjectStorage() {
    for (ChunkIterator b = chunks.begin(); b!=chunks.end(); ++b) {
      delete[] b->data;
    }
    chunks.clear();
}

#endif /* HUERTO_STORAGE_SMALL_OBJECT_STORAGE_HPP_ */

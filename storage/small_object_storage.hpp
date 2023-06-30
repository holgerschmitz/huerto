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
 * Removal of elements is performed by copying the last element of a chunk into the
 * position of the element to be removed and reducing the element count in the chunk by 
 * one.
 * 
 * Provides an STL style iterator and `begin()` and `end()` methods to iterate over all 
 * elements. Iteration is not in any guaranteed order.
 * 
 * @tparam T The type of values stored in the container
 * @tparam Traits A traits object that defines a `size_t chunkSize`
 */
template<typename T, class Traits = DefaultSmallObjectStorageTraits>
class SmallObjectStorage {
  public:
    /**
     * @brief The type of values stored in the container
     * 
     */
    typedef T value_type;
  private:
    /**
     * @brief Represents a fixed-size chunk of data that stores the objects
     * 
     * A chunk may not be fully filled with data either 
     * * if it is the last chunk in the list,
     * * or if an element has been removed after it was filled
     */
    struct DataChunk {
        /**
         * @brief Pointer to the array holding the elements.
         * 
         * The memory allocated for `data` is always given by the `chunkSize` trait.
         */
        T *data;

        /**
         * @brief The number of elements held by this chunk
         */
        size_t count;

        /**
         * @brief Do not default-construct a new Data Chunk object
         */
        DataChunk() = delete;

        /**
         * @brief Construct a new Data Chunk object given the data pointer 
         * 
         * @param data Pre-allocated pointer to the array holding the elements
         */
        DataChunk(T *data);

        /**
         * @brief Copy-construct a new Data Chunk object
         * 
         * A shallow copy will be made
         * 
         * @param chunk Another data chunk. 
         */
        DataChunk(const DataChunk &chunk);

        /**
         * @brief Add an element and return a reference to it
         * 
         * It must be guaranteed that the chunk is not full and can take the
         * additional element.
         * 
         * @return The element that was inserted
         */
        T &addElement();
    };

    /**
     * @brief An STL list of chunks
     */
    typedef std::list<DataChunk> ChunkList;

    /**
     * @brief An iterator over the STL list of chunks
     */
    typedef typename ChunkList::iterator ChunkIterator;

    /**
     * @brief A const iterator over the STL list of chunks
     */
    typedef typename ChunkList::const_iterator ChunkConstIterator;

    /**
     * @brief The STL list containing the chunks of data
     */
    ChunkList chunks;

    /**
     * @brief An iterator pointing to a chunk that may hold space for the next
     * element to be insterted 
     */
    ChunkIterator freeChunk;
  public:

    /**
     * @brief Construct a new Small Object Storage object
     */
    SmallObjectStorage() {
      freeChunk = chunks.end();
    }
    
    /**
     * @brief Do no copy-construct a new Small Object Storage object
     */
    SmallObjectStorage(const SmallObjectStorage&) = delete;

    /**
     * @brief Destroy the Small Object Storage object
     * 
     * Frees all the data associated with the chunks
     */
    ~SmallObjectStorage();

    /**
     * @brief A bidirectional iterator over the elements in the SmallObjectStorage
     */
    class iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
      private:
        friend class SmallObjectStorage<T, Traits>;
        /**
         * @brief The current position in the the chunk list
         */
        ChunkIterator chunkIter;

        /**
         * @brief The current position inside the chunk
         */
        size_t pos;

        /**
         * @brief Construct a new iterator object
         * 
         * @param chunkIter The current position in the the chunk list
         * @param pos The current position inside the chunk
         */
        iterator(ChunkIterator chunkIter, size_t pos = 0) :
          chunkIter(chunkIter), pos(pos) {}
      public:
        /**
         * @brief Construct a new iterator object
         */
        iterator() : pos(0) {}

        /**
         * @brief Copy-construct a new iterator object
         * 
         * @param it Another iterator
         */
        iterator(const iterator &it) = default;

        /**
         * @brief Prefix increment
         * 
         * @return The iterator
         */
        iterator& operator++() {
          if (++pos >= chunkIter->count) {
            ++chunkIter;
            pos = 0;
          }
          return *this;
        }

        /**
         * @brief Postfix increment
         * 
         * @return The iterator before the increment
         */
        iterator operator++(int) {
          iterator tmp(*this);
          operator++();
          return tmp;
        }

        /**
         * @brief Prefix decrement
         * 
         * @return The iterator
         */
        iterator& operator--() {
          if (pos == 0) {
            --chunkIter;
            pos = chunkIter->count - 1;
          } else {
            --pos;
          }
          return *this;
        }

        /**
         * @brief Postfix decrement
         * 
         * @return The iterator before the decrement
         */
        iterator operator--(int) {
          iterator tmp(*this);
          operator--();
          return tmp;
        }

        /**
         * @brief Comparison operator
         * 
         * @param rhs The iterator to compare to
         * @return True if both iterators point to the same element 
         */
        bool operator==(const iterator& rhs) const {
          return (chunkIter == rhs.chunkIter) && (pos == rhs.pos);
        }

        /**
         * @brief Inequality operator
         * 
         * @param rhs The iterator to compare to
         * @return True if both iterators point to the different elements
         */
        bool operator!=(const iterator& rhs) const {
          return (chunkIter != rhs.chunkIter) || (pos != rhs.pos);
        }

        /**
         * @brief Dereference operator
         * 
         * @return Reference to the current element
         */
        T& operator*() {
          return chunkIter->data[pos];
        }

        /**
         * @brief Dereference operator
         * 
         * @return Pointer to the current element
         */
        T* operator->() {
          return &(chunkIter->data[pos]);
        }

        /**
         * Assignment operator
         */
        iterator& operator=(const iterator& rhs) = default;
    };

    /**
     * @brief Returns an iterator to the beginning
     * 
     * @return An iterator to the beginning
     */
    iterator begin() {
      ChunkIterator b = chunks.begin();
      while (b->count == 0 && b != chunks.end()) {
        ++b;
      }
      return iterator(b);
    }

    /**
     * @brief Returns an iterator to the end
     * 
     * @return An iterator to the end
     */
    iterator end() {
      return iterator(chunks.end());
    }

    /**
     * @brief Add an element and return a reference to it
     * 
     * The reference is only valid until the next call to `removeElement`
     * 
     * @return Reference to the new element
     */
    T &addElement();

    /**
     * Remove an element from the storage.
     *
     * After deletion the iterator will point to the position after the deleted element.
     * This is in line with STL behaviour.
     *
     * Iterators pointing or references to elements before the deleted element remain valid.
     * Iterators pointing or references to elements after the deleted element become invalid.
     */
    iterator removeElement(const iterator&);

    /**
     * @brief Get the number of elements in the storage
     * 
     * @return the number of elements in the storage 
     */
    size_t getCount() const;

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
T &SmallObjectStorage<T, Traits>::DataChunk::addElement() {
  return data[count++];
}

template<typename T, class Traits>
T &SmallObjectStorage<T, Traits>::addElement() {
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
  return freeChunk->addElement();
}

template<typename T, class Traits>
typename SmallObjectStorage<T, Traits>::iterator SmallObjectStorage<T, Traits>::removeElement(const iterator &it_) {
  iterator it = it_;
  DataChunk &dchunk = *(it.chunkIter);
  T *data = dchunk.data;
  size_t pos = it.pos;
  size_t last = dchunk.count - 1;

  if (last>pos) std::swap(data[pos], data[last]);

  --(dchunk.count);
  // if no more elements in the chunk then delete the chunk
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
size_t SmallObjectStorage<T, Traits>::getCount() const
{
  size_t count = 0;
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

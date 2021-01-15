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

struct DefaultSmallObjectStorageTraits {
    static const size_t blockSize = 1000000;
};

template<typename T, class Traits = DefaultSmallObjectStorageTraits>
class SmallObjectStorage {
  public:
    typedef T value_type;
  private:
    struct DataBlock {
        T *data;
        size_t count;
        DataBlock() = delete;
        DataBlock(T *data);
        DataBlock(const DataBlock &block);
        T &addItem();
    };

    typedef std::list<DataBlock> BlockList;
    typedef typename BlockList::iterator BlockIterator;
    typedef typename BlockList::const_iterator BlockConstIterator;
    BlockList blocks;
    BlockIterator freeBlock;
  public:

    SmallObjectStorage() {
//      blocks.push_front(DataBlock(new T[Traits::blockSize]));
      freeBlock = blocks.end();
    }

    ~SmallObjectStorage();

    class iterator : public std::iterator<std::bidirectional_iterator_tag, T> {
      private:
        friend class SmallObjectStorage<T, Traits>;
        BlockIterator blockIter;
        size_t pos;

        iterator(BlockIterator blockIter_, long pos_ = 0) :
          blockIter(blockIter_), pos(pos_) {}
      public:
        iterator() : pos(0) {}

        iterator(const iterator &it) : blockIter(it.blockIter), pos(it.pos) {}

        iterator& operator++() {
          if (++pos >= blockIter->count) {
            ++blockIter;
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
            --blockIter;
            pos = blockIter->count - 1;
          }
          return *this;
        }

        iterator operator--(int) {
          iterator tmp(*this);
          operator--();
          return tmp;
        }

        bool operator==(const iterator& rhs) const {
          return (blockIter == rhs.blockIter) && (pos == rhs.pos);
        }

        bool operator!=(const iterator& rhs) const {
          return (blockIter != rhs.blockIter) || (pos != rhs.pos);
        }

        T& operator*() {
          return blockIter->data[pos];
        }

        T* operator->() {
          return &(blockIter->data[pos]);
        }
    };

    iterator begin() {
      BlockIterator b = blocks.begin();
      while (b->count == 0 && b != blocks.end()) {
        ++b;
      }
      return iterator(b);
    }

    iterator end() {
      return iterator(blocks.end());
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

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::DataBlock::DataBlock(T *data) : data(data), count(0) {
//    std::cerr << "Created datablock " << data << std::endl;
}

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::DataBlock::DataBlock(const DataBlock &block) : data(block.data), count(block.count) {
//    std::cerr << "Copied datablock " << data << std::endl;
}

template<typename T, class Traits>
T &SmallObjectStorage<T, Traits>::DataBlock::addItem() {
  return data[count++];
}

template<typename T, class Traits>
T &SmallObjectStorage<T, Traits>::addItem() {
  // find a free block or create one
  if ((freeBlock == blocks.end()) || (Traits::blockSize == freeBlock->count)) {
    freeBlock = blocks.begin();

    while((freeBlock != blocks.end()) && (Traits::blockSize <= freeBlock->count)) {
      ++freeBlock;
    }

    if (freeBlock == blocks.end()) {
//      std::cerr << "New DataBlock" << std::endl;
      freeBlock = blocks.insert(blocks.begin(), DataBlock(new T[Traits::blockSize]));
    }
  }
  return freeBlock->addItem();
}

template<typename T, class Traits>
typename SmallObjectStorage<T, Traits>::iterator SmallObjectStorage<T, Traits>::removeItem(const iterator &it_) {
  iterator it = it_;
  DataBlock &dblock = *(it.blockIter);
  T *data = dblock.data;
  long pos = it.pos;
  long last = dblock.count - 1;

  if (last>pos) std::swap(data[pos], data[last]);

  --(dblock.count);
  // if no more rays in the block then delete the block
  if (0 == dblock.count && blocks.size() > 1) {
    std::cerr << "Deleting data" << dblock.data << std::endl;
    delete[] dblock.data;
    it.blockIter = blocks.erase(it.blockIter);
    it.pos = 0;
    freeBlock = blocks.begin();
  } else if (last==pos) {
    ++it.blockIter;
    it.pos = 0;
  }

  // for performance, see if there are more free spaces in this block than in the
  // freeBlock, and update freeBlock accordingly
//  if ((it.blockIter != blocks.end()) && (it.blockIter->count < freeBlock->count))
//  {
//    freeBlock = it.blockIter;
//  }
  return it;
}

template<typename T, class Traits>
long SmallObjectStorage<T, Traits>::getCount() const
{
  long count = 0;
  for (BlockConstIterator b = blocks.begin(); b!=blocks.end(); ++b)
    count += b->count;
  return count;
}

template<typename T, class Traits>
void SmallObjectStorage<T, Traits>::clear() {
    for (BlockIterator b = blocks.begin(); b!=blocks.end(); ++b) {
      delete[] b->data;
    }
    blocks.clear();
    freeBlock = blocks.end();
}

template<typename T, class Traits>
SmallObjectStorage<T, Traits>::~SmallObjectStorage() {
    for (BlockIterator b = blocks.begin(); b!=blocks.end(); ++b) {
      delete[] b->data;
    }
    blocks.clear();
}
#endif /* HUERTO_STORAGE_SMALL_OBJECT_STORAGE_HPP_ */

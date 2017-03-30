/* 
 *  Copyright (c) 2011 Daisuke Okanohara
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */


#ifndef RANK_VECTOR_HPP_
#define RANK_VECTOR_HPP_

#include <vector>
#include <stdint.h>

namespace dag{

/**
 * Bit Vector supporing Rank operation
 */
class rank_vector{
public:

  /**
   * Constructor
   */
  rank_vector(): size_(0), one_num_(0){
    bits_.push_back(0);
    lblocks_.push_back(0);
    sblocks_.push_back(0);
  }

  /**
   * Destructor
   */
  ~rank_vector(){
  }

  /**
   * Add bit to the end of the vector
   * @param bit a bit to be added
   */
  void push_back(uint64_t bit){
    if (bit){
      bits_[size_ / BLOCKSIZE] |= (1LLU << (size_ % BLOCKSIZE)); 
    }
    size_++;
    if ((size_ % BLOCKSIZE) == 0){
      add_block();
    }
  }

  /**
   * Get the pos-th bit
   * @param pos the index
   * @return the pos-th bit
   */
  uint64_t get_bit(uint64_t pos) const{
    return (bits_[pos/BLOCKSIZE] >> (pos % BLOCKSIZE)) & 0x1LLU; 
  }

  /**
   * Calculate the number of ones in bits_[0...pos-1] in O(1) time.
   * @param pos the position in the bit array
   * @return the number of ones in bits_[0...pos-1]
   */
  uint64_t rank(uint64_t pos) const{
    return lblocks_[pos/LBLOCKSIZE] 
      + sblocks_[pos/BLOCKSIZE] 
      + pop_count(bits_[pos/BLOCKSIZE] & ((1LLU << (pos % BLOCKSIZE)) - 1));
   }

  /**
   * Return the size of bit array in bits.
   * @return the number of bits
   */
  uint64_t size() const{
    return size_;
  }

  /**
   * Swap the content in bit vector
   * @param rv the rank_vector to be swapped
   */
  void swap(rank_vector& rv){
    bits_.swap(rv.bits_);
    lblocks_.swap(rv.lblocks_);
    sblocks_.swap(rv.sblocks_); 
    std::swap(size_, rv.size_);
    std::swap(one_num_, rv.one_num_);
  }

 private:
  static const uint64_t LBLOCKSIZE = 256;
  static const uint64_t BLOCKSIZE = 64;

  void add_block(){
    if (bits_.size() > 0){
      one_num_ += pop_count(bits_.back());
    }

    if (size_ % LBLOCKSIZE == 0){
      lblocks_.push_back(one_num_);
    }
    sblocks_.push_back(one_num_ - lblocks_[size_ / LBLOCKSIZE]);
    bits_.push_back(0LLU);
  }

  inline static uint64_t pop_count(uint64_t x){
    x = x - ((x & 0xAAAAAAAAAAAAAAAALLU) >> 1);
    x = (x & 0x3333333333333333LLU) + ((x >> 2) & 0x3333333333333333LLU);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
    return x * 0x0101010101010101LLU >> 56;
  }

  std::vector<uint64_t> bits_;    /// bit array
  std::vector<uint64_t> lblocks_; /// rank results for large blocks
  std::vector<uint8_t> sblocks_;  /// rank results for small blocks
  uint64_t size_;                 /// the length of bit array
  uint64_t one_num_;              /// the number of ones in the bit array
};


}

#endif // RANK_VECTOR_HPP_

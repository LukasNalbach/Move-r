/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 // * @file StepCode.hpp
 // * @brief Simple space economic code for uints in which each uint is stored in 4*i bits with smallest possible integer i.
 // * @author Tomohiro I
 // * @date 2017-01-24
 */
#ifndef INCLUDE_GUARD_StepCode
#define INCLUDE_GUARD_StepCode

#include <stdint.h> // include uint64_t etc.
#include <assert.h>

#include <iostream>
#include <iterator>
#include <algorithm>

#include "BitsUtil.hpp"
#include "MemUtil.hpp"

namespace itmmti
{
  /*!
   * @brief Utilities for StepCode.
   * @par notation
   *   Each uint is stored in 4*i bits with smallest possible integer i.
   *   Suppose that i-th (0base) uint occupies w_i bits (remark that w_i in {4, 8, 12, 16, ..., 64}).
   *   w_i/4 - 1 is called wCode of w_i (remark that w_i/4 - 1 in {0, 1, 2, 3, ..., 7}). wCodes are stored in wCodes_ using 4 bits each.
   */
  struct StepCodeUtil
  {
    static constexpr uint8_t kStep{4}; // Bits of each step.
    static constexpr uint8_t kWCBits{4}; // Bits of each wCode.
    static constexpr uint8_t kWCNum{16}; // Number of wCodes in 64 bits.


    /*!
     * @brief Calculate wCode for a given value.
     */
    static uint8_t calcWCode
    (
     const uint64_t val //!< in [0, 2^{64}).
     ) noexcept {
      return (bits::bitSize(val) - 1) / kStep;
    }


    /*!
     * @brief Calculate smallest bit-width in {4, 8, 12, ..., 64} to store a given value.
     */
    static uint8_t calcSteppedW
    (
     const uint64_t val //!< in [0, 2^{64}).
     ) noexcept {
      return (bits::bitSize(val) + 3) / kStep * kStep;
    }


    /*!
     * @brief Calculate wCode from steppedW.
     */
    static uint8_t calcWCodeFromSteppedW
    (
     const uint8_t steppedW //!< in {4, 8, 12, ..., 64}.
     ) noexcept {
      assert(steppedW <= 64);
      assert(steppedW % 4 == 0);

      return (steppedW / kStep) - 1;
    }


    /*!
     * @brief Read bit-width of "idx"-th (0base) values.
     */
    static uint8_t readW
    (
     const uint64_t * wCodes, //!< Array storing wCodes using 4 bits each.
     const uint64_t idx //!< "(idx + 1) * 4 - 1"-th bit of wCodes array should not be out-of-bound.
     ) noexcept {
      return kStep * static_cast<uint8_t>(bits::readWBits_S(wCodes, kWCBits * idx, ctcbits::UINTW_MAX(kWCBits)) + 1);
    }


    /*!
     * @brief Write wCode for "idx"-th (0base) value.
     */
    static void writeWCode
    (
     const uint8_t wCode, //!< in [0..2^4).
     uint64_t * wCodes, //!< Array to which "wCode" is written (written at bit-region "[[idx * 4..(idx + 1) * 4))").
     const uint64_t idx //!< "(idx + 1) * 4 - 1"-th bit of wCodes array should not be out-of-bound.
     ) noexcept {
      assert(wCode <= ctcbits::UINTW_MAX(kWCBits));

      bits::writeWBits_S(wCode, wCodes, kWCBits * idx, ctcbits::UINTW_MAX(kWCBits));
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value.
     */
    static uint64_t calcBitPos
    (
     const uint64_t * wCodes, //!< Array storing wCodes using 4 bits each.
     const uint64_t idx //!< "(idx + 1) * 4 - 1"-th bit of wCodes array should not be out-of-bound.
     ) noexcept {
      auto sum = idx;
      for (uint64_t i = 0; i < idx / kWCNum; ++i) {
        sum += StepCodeUtil::sumWCodes(wCodes[i]);
      }
      uint8_t rem = idx % kWCNum;
      if (rem) {
        sum += StepCodeUtil::sumWCodes(wCodes[idx / kWCNum] & bits::UINTW_MAX(rem * kWCBits));
      }
      return sum * kStep;
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value with auxiliary array of middle level.
     */
    static uint64_t calcBitPos
    (
     const uint64_t * wCodes, //!< Array storing wCodes using 4 bits each.
     const uint64_t idx, //!< "(idx + 1) * 4 - 1"-th bit of wCodesArray should not be out-of-bound.
     const uint8_t * wCodesAuxM //! Auxiliary array: "wCodesAuxM[i]" stores sum of wCodes in "wCodes[i]".
     ) noexcept {
      auto sum = idx;
      for (uint64_t i = 0; i < idx / kWCNum; ++i) {
        sum += wCodesAuxM[i];
      }
      uint8_t rem = idx % kWCNum;
      if (rem) {
        sum += StepCodeUtil::sumWCodes(wCodes[idx / kWCNum] & bits::UINTW_MAX(rem * kWCBits));
      }
      return sum * kStep;
    }


    /*!
     * @brief Calculate sum of bit-widths of values indexed from "beg" (inclusive) to "end" (exclusive).
     * @pre "beg <= end".
     */
    static uint64_t sumW
    (
     const uint64_t * wCodes, //!< Array storing wCodes using 4 bits each.
     const uint64_t beg, //!< beg (<= end).
     const uint64_t end //!< "(end + 1) * 4 - 1"-th bit of wCodes array should not be out-of-bound.
     ) noexcept {
      assert(beg <= end);

      uint64_t sum = 0;
      for (uint64_t i = beg; i < end; ++i) {
        sum += readW(wCodes, i);
      }
      return sum;
    }


    /*!
     * @brief Calculate sum of wCodes in a given 64bits uint.
     */
    static uint8_t sumWCodes
    (
     uint64_t wCodesInWord //!< 64bits uint storing wCodes.
     ) noexcept {
      wCodesInWord = ((wCodesInWord & 0xf0f0f0f0f0f0f0f0ULL)>>4) + (wCodesInWord & 0x0f0f0f0f0f0f0f0fULL);
      wCodesInWord += (wCodesInWord>>8);
      wCodesInWord += (wCodesInWord>>16);
      wCodesInWord += (wCodesInWord>>32);
      return wCodesInWord & 0xff;
    }


    /*!
     * @brief Change wCode array by inserting and/or deleting wCodes.
     * @pre The resulting size should be within capacity.
     * @node
     *   This function does not touch bit array storing values.
     *   If there is a need to shift bit positions of values, call changeValPos().
     *   If there is a need to update wCodesAuxM array, call updateWCodesAuxM().
     */
    // template <typename SizeT>
    // static void changeWCodes
    // (
    //  const uint64_t * srcWCodes, //!< Array storing source of wCodes using 4 bits each.
    //  const uint64_t srcIdxBeg, //!< Beginning idx of src.
    //  const uint64_t srcLen, //!< Length of wCodes of src to insert.
    //  uint64_t * tgtWCodes, //!< Target wCodes array to change.
    //  const uint64_t tgtIdxBeg, //!< Beginning idx of tgt.
    //  const uint64_t tgtLen, //!< Length of wCodes of tgt to delete.
    //  SizeT & size
    //  ) noexcept {
    //   assert(tgtIdxBeg + tgtLen <= size);
    //   assert(size + srcLen >= tgtLen);

    //   const auto tailNum = size - (tgtIdxBeg + tgtLen); // at least 0 by assumption.
    //   if (srcLen != tgtLen && tailNum) { // Need to shift elements in wCodes_.
    //     const auto mvSrcPos = (tgtIdxBeg + tgtLen) * StepCodeUtil::kWCBits;
    //     const auto mvTgtPos = (tgtIdxBeg + srcLen) * StepCodeUtil::kWCBits;
    //     bits::mvBits(tgtWCodes + (mvSrcPos / 64), mvSrcPos % 64, tgtWCodes + (mvTgtPos / 64), mvTgtPos % 64, tailNum * StepCodeUtil::kWCBits);
    //   }
    //   if (srcLen) {
    //     const auto mvSrcPos = srcIdxBeg * StepCodeUtil::kWCBits;
    //     const auto mvTgtPos = tgtIdxBeg * StepCodeUtil::kWCBits;
    //     bits::mvBits(srcWCodes + (mvSrcPos / 64), mvSrcPos % 64, tgtWCodes + (mvTgtPos / 64), mvTgtPos % 64, srcLen * StepCodeUtil::kWCBits);
    //   }
    //   size += srcLen - tgtLen;
    // }


    /*!
     * @brief Update wCodesAuxM.
     */
    static void updateWCodesAuxM
    (
     const uint64_t * wCodes, //!< Array storing wCodes using 4 bits each.
     uint8_t * wCodesAuxM, //!< Auxiliary array: "wCodesAuxM[i]" stores sum of wCodes in "wCodes[i]".
     const uint64_t beg_word, //!< beg index of 64bits uint array of wCodes
     const uint64_t end_word //!< end (exclusive) index of 64bits uint array of wCodes
     ) noexcept {
      assert(beg_word <= end_word);

      for (uint64_t i = beg_word; i < end_word; ++i) {
        wCodesAuxM[i] = StepCodeUtil::sumWCodes(wCodes[i]);
      }
    }


    /*!
     * @brief Change (shift) value positions.
     * @pre The resulting bit size should be within capacity.
     */
    // template <typename SizeT>
    // static void changeValPos
    // (
    //  uint64_t * vals, //!< Array storing values.
    //  SizeT & bitSize, //!< [in,out] bit size.
    //  const uint64_t bitPos, //!< Beginning bit-pos in vals.
    //  const uint64_t insBitLen, //!< Bit-length to insert.
    //  const uint64_t delBitLen //!< Bit-length to delete.
    //  ) noexcept {
    //   const uint64_t srcPos = bitPos + delBitLen;
    //   const uint64_t tgtPos = bitPos + insBitLen;
    //   bits::mvBits(vals + (srcPos / 64), srcPos % 64, vals + (tgtPos / 64), tgtPos % 64, bitSize - srcPos);
    //   bitSize += insBitLen - delBitLen;
    // }


    /*!
     * @brief Move vals.
     * @pre Bit-regions should be within capacity.
     */
    static void mvVals
    (
     const uint64_t * srcVals, //!< Array storing source values.
     const uint64_t srcBitPos, //!< Beginning bit-pos in "srcVals".
     uint64_t * tgtVals, //!< Target array.
     const uint64_t tgtBitPos, //!< Beginning bit-pos in "tgtVals".
     const uint64_t bitLen //!< Bit-len to move.
     ) noexcept {
      bits::mvBits(srcVals, srcBitPos, tgtVals, tgtBitPos, bitLen);
    }
  };

  




  /////////////////////////////////////////////////////////////////////////
  /*!
   * @brief Core part of StepCode.
   * @tparam
   *   kMaxNum: Max number of uints to be stored. It is fixed as StepCode is intended to be used for small kMaxNum.
   * @attention
   *   We do not provide a constant-time random access.
   *   bitSize and bitCapacity of "vals_" should be maintained manually.
   * @par notation
   *   Suppose that i-th (0base) uint occupies w_i bits (remark that w_i in {4, 8, 12, 16, ..., 64}).
   *   w_i/4 - 1 is called wCode of w_i (remark that w_i/4 - 1 in {0, 1, 2, 3, ..., 7}). wCodes are stored in wCodes_ using 4 bits each.
   */
  template<uint32_t kMaxNum>
  class StepCodeCore
  {
  private:
    uint64_t wCodes_[kMaxNum / StepCodeUtil::kWCNum]; //!< Store wCodes using 4 bits each.
    uint64_t * vals_; //!< Array to store values.


  public:
    StepCodeCore() : vals_(nullptr)
    {
    }


    template<typename SizeT>
    StepCodeCore
    (
     SizeT & bitCapacity, //!< [in,out] Give bit capacity. It will be modified.
     const size_t initBitCapacity = 0 //!< Initial capacity of bv_.
     ) : vals_(nullptr) {
      assert(initBitCapacity <= ctcbits::UINTW_MAX(58));

      changeValsBitCapacity(bitCapacity, 0, initBitCapacity);
    }


    StepCodeCore(const StepCodeCore & other) = delete;
    StepCodeCore & operator= (const StepCodeCore & other) = delete;
    StepCodeCore(const StepCodeCore && other) = delete;
    StepCodeCore & operator= (const StepCodeCore && other) = delete;


    /*!
     * @brief Constructor for copy.
     * @attention
     *   Since the contents of "other" are copied, it may take time when other.size_ is large.
     */
    template<typename SizeT>
    StepCodeCore
    (
     const StepCodeCore & other, //!< Other StepCodeCore object.
     const SizeT otherBitSize, //!< Bit size of "other".
     const SizeT otherSize, //!< Size of "other".
     SizeT & bitCapacity, //!< [out] Capture valBitCapacity.
     SizeT & bitSize, //!< [out] Capture valBitSize.
     SizeT & size //!< [out] Capture size.
     ) : vals_(nullptr) {
      this->cpStepCodeCore(other, otherBitSize, otherSize, bitCapacity, bitSize, size);
    }


    template<typename SizeT>
    void cpStepCodeCore
    (
     const StepCodeCore & other, //!< Other StepCodeCore object.
     const SizeT otherBitSize, //!< Bit size of "other".
     const SizeT otherSize, //!< Size of "other".
     SizeT & bitCapacity, //!< [out] Capture valBitCapacity.
     SizeT & bitSize, //!< [out] Capture valBitSize.
     SizeT & size //!< [out] Capture size.
     ) {
      if (otherSize > 0) { // Lazy reservation: If size_ == 0, we do not reserve anything.
        bitCapacity = static_cast<SizeT>(this->setBitCapacity(static_cast<size_t>(std::max(otherBitSize, bitSize))));
        bits::cpBytes(other.vals_, this->vals_, (otherBitSize + 7) / 8);
        bits::cpBytes(other.wCodes_, this->wCodes_, (otherSize * StepCodeUtil::kWCBits + 7) / 8);
        size = otherSize;
        bitSize = otherBitSize;
      } else {
        size = bitSize = 0;
      }
    }


    //// Move
    /*!
     * @brief Constructor for move.
     * @attention
     *   "other" is initialized to an object with capacity = 0.
     */
    template<typename SizeT>
    StepCodeCore
    (
     StepCodeCore && other, //!< Other StepCodeCore object.
     const SizeT otherBitCapacity, //!< Bit capacity of "other".
     const SizeT otherBitSize, //!< Bit size of "other".
     const SizeT otherSize, //!< Size of "other".
     SizeT & bitCapacity, //!< [out] Capture valBitCapacity.
     SizeT & bitSize, //!< [out] Capture valBitSize.
     SizeT & size //!< [out] Capture size.
     ) : vals_(other.vals_) {
      this->mvStepCodeCore(other, otherBitCapacity, otherBitSize, otherSize, bitCapacity, bitSize, size);
    }


    template<typename SizeT>
    void mvStepCodeCore
    (
     StepCodeCore && other, //!< Other StepCodeCore object.
     const SizeT otherBitCapacity, //!< Bit capacity of "other".
     const SizeT otherBitSize, //!< Bit size of "other".
     const SizeT otherSize, //!< Size of "other".
     SizeT & bitCapacity, //!< [out] Capture valBitCapacity.
     SizeT & bitSize, //!< [out] Capture valBitSize.
     SizeT & size //!< [out] Capture size.
     ) {
      vals_ = other.vals_;
      bitCapacity = otherBitCapacity;
      bitSize = otherBitSize;
      size = otherSize;
      if (size) {
        bits::cpBytes(other.wCodes_, wCodes_, (size * StepCodeUtil::kWCBits + 7) / 8);
      }
      other.vals_ = nullptr;
    }


    ~StepCodeCore()
    {
      free(vals_);
    }


    /*!
     * @brief Get read-only vals_ array pointer.
     */
    const uint64_t * getConstPtr_vals() const noexcept
    {
      return vals_;
    }


    /*!
     * @brief Get read-only wCodes_ array pointer.
     */
    const uint64_t * getConstPtr_wCodes() const noexcept
    {
      return wCodes_;
    }


    /*!
     * @brief Calculate sum of bit-widths of values indexed from "beg" (inclusive) to "end" (exclusive).
     * @pre "beg <= end (<= size)". Note that "size" is maintained outside of StepCodeCore.
     */
    uint64_t sumW
    (
     const uint64_t beg, //!< in [0..size].
     const uint64_t end //!< in [beg..size].
     ) const noexcept {
      assert(beg <= end);

      return StepCodeUtil::sumW(wCodes_, beg, end);
    }


    /*!
     * @brief Read bit-width of "idx"-th (0base) values.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     */
    uint8_t readW
    (
     const uint64_t idx //!< in [0..size).
     ) const noexcept {
      return StepCodeUtil::readW(wCodes_, idx);
    }


    /*!
     * @brief Write wCode for "idx"-th (0base) value.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     */
    void writeWCode
    (
     const uint8_t wCode, //!< in [0..2^4).
     const uint64_t idx //!< in [0..size).
     ) noexcept {
      assert(wCode <= ctcbits::UINTW_MAX(StepCodeUtil::kWCBits));

      StepCodeUtil::writeWCode(wCode, wCodes_, idx);
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in bv.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     */
    uint64_t calcBitPos
    (
     const uint64_t idx //!< in [0..size).
     ) const noexcept {
      return StepCodeUtil::calcBitPos(wCodes_, idx);
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in bv.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     */
    uint64_t calcBitPos
    (
     const uint64_t idx, //!< in [0..size).
     const uint8_t * wCodesAuxM //! Auxiliary array: "wCodesAuxM[i]" stores sum of wCodes in "wCodes[i]".
     ) const noexcept {
      return StepCodeUtil::calcBitPos(wCodes_, idx, wCodesAuxM);
    }


    /*!
     * @brief Read "w"-bits written in the bit-region beginning at array_[[bitPos..]] in bv_.
     * @return Value represented by array_[[bitPos..bitPos+w)).
     * @pre The bit-region must not be out of bounds.
     */
    uint64_t readWBits
    (
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of bit-region.
     const uint8_t steppedW //!< Bit-width in {4, 8, 12, 16, ..., 64}.
     ) const noexcept {
      assert(steppedW <= 64);
      assert(steppedW % 4 == 0);

      return bits::readWBits(vals_, bitPos, steppedW, bits::UINTW_MAX(steppedW));
    }


    /*!
     * @brief Write "w"-bit value "val" to the bit-region beginning at array_[[bitPos..]].
     * @pre The bit-region must not be out of bounds.
     */
    void writeWBits
    (
     const uint64_t val, //!< in [0..2^steppedW).
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of bit-region.
     const uint8_t steppedW //!< Bit-width in {4, 8, 12, 16, ..., 64}.
     ) noexcept {
      assert(steppedW <= 64);
      assert(steppedW % 4 == 0);
      assert(bits::bitSize(val) <= steppedW);

      bits::writeWBits(val, vals_, bitPos, steppedW, bits::UINTW_MAX(steppedW));
    }


    /*!
     * @brief Read "idx"-th (0base) value.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     * @note It may be inefficient for sequencial access. Use readW and readWBits in that case.
     */
    uint64_t read
    (
     uint64_t idx //!< in [0, size_).
     ) const noexcept {
      const auto bitPos = StepCodeUtil::calcBitPos(wCodes_, idx);
      const auto w = readW(idx);
      return readWBits(bitPos, w);
    }


    /*!
     * @brief Read "idx"-th (0base) value.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     * @note It may be inefficient for sequencial access. Use readW and readWBits in that case.
     */
    uint64_t read
    (
     uint64_t idx, //!< in [0, size_).
     const uint8_t * wCodesAuxM //! Auxiliary array: "wCodesAuxM[i]" stores sum of wCodes in "wCodes[i]".
     ) const noexcept {
      const auto bitPos = StepCodeUtil::calcBitPos(wCodes_, idx, wCodesAuxM);
      const auto w = readW(idx);
      return readWBits(bitPos, w);
    }


    /*!
     * @brief Read "idx"-th (0base) value.
     * @note It may be inefficient for sequencial access. Use readW and readWBits in that case.
     */
    uint64_t read_naive
    (
     uint64_t idx //!< in [0, size_).
     ) const noexcept {
      uint64_t bitPos = 0;
      for (uint64_t i = 0; i < idx; ++i) {
        auto bitPos = readW(i);
      }
      const auto w = readW(idx);
      return readWBits(bitPos, w);
    }


    /*!
     * @brief Write uint "val" at the end.
     * @pre The resulting size_ and bitSize_ should be within capacity.
     */
    // template <typename SizeT>
    // void append
    // (
    //  const uint64_t val, //!< Uint val.
    //  const uint8_t steppedW, //!<  bit-width in {4, 8, 12, ..., 64} into which val fits.
    //  SizeT & size, //!< [in,out] Size.
    //  SizeT & bitSize //!< [in,out] Bit size.
    //  ) {
    //   assert(size < kMaxNum);
    //   assert(steppedW <= 64 && steppedW % 4 == 0);
    //   assert(bits::bitSize(val) <= steppedW);

    //   const uint8_t wCode = steppedW / StepCodeUtil::kStep - 1;
    //   const auto begPos = bitSize;
    //   bitSize += steppedW;
    //   writeWBits(val, bitSize, steppedW);
    //   bits::writeWBits_S(wCode, wCodes_, StepCodeUtil::kWCBits * size, ctcbits::UINTW_MAX(StepCodeUtil::kWCBits));
    //   ++size;
    // }


    /*!
     * @brief Write uint "val" at the end.
     * @pre The resulting size_ and bitSize_ should be within capacity.
     */
    // template <typename SizeT>
    // void append
    // (
    //  const uint64_t val, //!< Uint val.
    //  SizeT & size, //!< [in,out] Size.
    //  SizeT & bitSize //!< [in,out] Bit size.
    //  ) {
    //   assert(size < kMaxNum);

    //   const auto wCode = StepCodeUtil::calcWCode(val);
    //   const auto w = (wCode + 1) * StepCodeUtil::kStep;
    //   const auto begPos = bitSize;
    //   bitSize += w;
    //   writeWBits(val, begPos, w);
    //   bits::writeWBits_S(wCode, wCodes_, StepCodeUtil::kWCBits * size, ctcbits::UINTW_MAX(StepCodeUtil::kWCBits));
    //   ++size;
    // }


    /*!
     * @brief Rewrite value at "idx" to "val" without changing bit-width.
     * @pre "idx < size". Note that "size" is maintained outside of StepCodeCore.
     */
    void rewriteVal
    (
     const uint64_t val, //!< should fit in readW(idx) bits.
     const uint64_t idx, //!< in [0, size).
     const uint64_t bitPos //!< should be bit-position of "idx"-th (0base) values in bv_.
     ) {
      assert(bits::bitSize(val) <= UINTW_MAX(readW(idx)));

      const auto w = readW(idx);
      writeWBits(val, bitPos, w);
    }


    /*!
     * @brief Move wCodes to "wCodes_"
     * @pre Be sure to work within capacity of "srcWCodes" and "wCodes_".
     * @node
     *   - Be sure to maintain size outside if needed.
     *   - Be sure maintain value positions outside if needed.
     */
    void mvWCodes
    (
     const uint64_t * srcWCodes, //!< Array storing source of wCodes using 4 bits each.
     const uint64_t srcIdxBeg, //!< Beginning idx of src.
     const uint64_t tgtIdxBeg, //!< Beginning idx of tgt.
     const uint64_t len //!< Length of wCodes to move.
     ) noexcept {
      const auto mvSrcPos = srcIdxBeg * StepCodeUtil::kWCBits;
      const auto mvTgtPos = tgtIdxBeg * StepCodeUtil::kWCBits;
      bits::mvBits(srcWCodes, mvSrcPos, wCodes_, mvTgtPos, len * StepCodeUtil::kWCBits);
    }


    /*!
     * @brief Move vals.
     * @pre Bit-regions should be within capacity.
     * @node
     *   - Be sure to maintain bitSize outside if needed.
     */
    void mvVals
    (
     const uint64_t * srcVals, //!< Array storing source values.
     const uint64_t srcBitPos, //!< Beginning bit-pos in "srcVals".
     const uint64_t tgtBitPos, //!< Beginning bit-pos in "tgtVals".
     const uint64_t bitLen //!< Bit-len to move.
     ) noexcept {
      bits::mvBits(srcVals, srcBitPos, vals_, tgtBitPos, bitLen);
    }


    /*!
     * @brief Change wCode array by inserting and deleting wCodes.
     * @pre The resulting size should be within capacity.
     * @node
     *   If (tailNum && insBitLen != delBitLen), changeValPos should be called to maintain (shift) value positions.
     */
    // template <typename SizeT>
    // void changeWCodes
    // (
    //  const uint64_t * srcWCodes, //!< Array storing source of wCodes using 4 bits each.
    //  const uint64_t srcIdxBeg, //!< Beginning idx of src.
    //  const uint64_t srcLen, //!< Length of wCodes of src to insert.
    //  const uint64_t tgtIdxBeg, //!< Beginning idx of tgt.
    //  const uint64_t tgtLen, //!< Length of wCodes of tgt to delete.
    //  SizeT & size
    //  ) noexcept {
    //   assert(tgtIdxBeg + tgtLen <= size);
    //   assert(size + srcLen >= tgtLen);

    //   StepCodeUtil::changeWCodes(srcWCodes, srcIdxBeg, srcLen, wCodes_, tgtIdxBeg, tgtLen, size);
    // }


    /*!
     * @brief Update wCodesAuxM.
     */
    void updateWCodesAuxM
    (
     uint8_t * wCodesAuxM, //!< Auxiliary array: "wCodesAuxM[i]" stores sum of wCodes in "wCodes[i]".
     const uint64_t beg_word, //!< beg index of 64bits uint array of wCodes
     const uint64_t end_word //!< end (exclusive) index of 64bits uint array of wCodes
     ) noexcept {
      for (uint64_t i = beg_word; i < end_word; ++i) {
        wCodesAuxM[i] = StepCodeUtil::sumWCodes(wCodes_[i]);
      }
    }


    /*!
     * @brief Change (shift) value positions.
     * @pre The resulting size should be within capacity.
     */
    // template <typename SizeT>
    // void changeValPos
    // (
    //  SizeT & bitSize,
    //  const uint64_t bitPos, //!< bitPos of bv_.
    //  const uint64_t insBitLen, //!< Bit-length to insert in bv_
    //  const uint64_t delBitLen //!< Bit-length to delete in bv_
    //  ) noexcept {
    //   const uint64_t srcPos = bitPos + delBitLen;
    //   const uint64_t tgtPos = bitPos + insBitLen;
    //   bits::mvBits(vals_ + (srcPos / 64), srcPos % 64, vals_ + (tgtPos / 64), tgtPos % 64, bitSize - srcPos);
    //   bitSize += insBitLen - delBitLen;
    // }


    /*!
     * @brief Get capacity (fixed).
     */
    size_t capacity() const noexcept {
      return kMaxNum;
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     * @node It does not care size.
     */
    size_t setBitCapacity
    (
     const size_t givenCapacity
     ) {
      assert(givenCapacity <= ctcbits::UINTW_MAX(58));

      const size_t newLen = (givenCapacity + 63) / 64; // +63 for roundup
      if (newLen > 0) {
        memutil::realloc_AbortOnFail<uint64_t>(vals_, newLen);
        return newLen * 64;
      } else {
        memutil::safefree(vals_);
        return 0;
      }
    }
  };






  /////////////////////////////////////////////////////////////////////////
  /*!
   * @brief Simple space economic code for uints in which each uint is stored in 4*i bits with smallest possible integer i.
   * @tparam
   *   kMaxNum: Max number of uints to be stored. It is fixed as StepCode is intended to be used for small kMaxNum.
   * @attention
   *   We do not provide a constant-time random access.
   * @par notation
   *   Suppose that i-th (0base) uint occupies w_i bits (remark that w_i in {4, 8, 12, 16, ..., 64}).
   *   w_i/4 - 1 is called wCode of w_i (remark that w_i/4 - 1 in {0, 1, 2, 3, ..., 7}). wCodes are stored in wCodes_ using 4 bits each.
   */
  template<uint32_t kMaxNum, typename SizeT = uint32_t>
  class StepCode
  {
  private:
    //// Private member variables.
    StepCodeCore<kMaxNum> core_;
    SizeT bitCapacity_; //!< Current capacity (must be in [0, 2^58)).
    SizeT bitSize_; //!< Current size (must be in [0, capacity_]).
    SizeT size_; //!< Current size (number of elements).


  public:
    StepCode
    (
     size_t initBitCapacity = 0 //!< Initial capacity of bv_.
     ) : core_(), bitCapacity_(0), bitSize_(0), size_(0) {
      assert(initBitCapacity <= ctcbits::UINTW_MAX(sizeof(SizeT) * 8));
      assert(initBitCapacity <= ctcbits::UINTW_MAX(58));

      changeBitCapacity(initBitCapacity);
    }


    //// Copy
    /*!
     * @brief Copy constructor.
     * @attention
     *   Since the contents of "other" are copied, it may take time when other.size_ is large.
     */
    StepCode
    (
     const StepCode & other
     ) : core_(), bitCapacity_(0), bitSize_(other.bitSize_), size_(other.size_) {
      core_.cpStepCodeCore(other.core_, other.bitSize_, other.size_, bitCapacity_, bitSize_, size_);
    }


    /*!
     * @brief Assignment operator.
     * @attention
     *   - The contents of 'other' are copied.
     *   - It may take time when other.size_ is large.
     *   - The original contents are discarded.
     *   - Lazy allocation: Changing capacity is done only if it is needed.
     */
    StepCode & operator=
    (
     const StepCode & other
     ) {
      if (this != &other) {
        core_.cpStepCodeCore(other.core_, other.bitSize_, other.size_, bitCapacity_, bitSize_, size_);
      }
      return *this;
    }


    //// Move
    /*!
     * @brief Move constructor.
     * @attention
     *   "other" is initialized to an object with capacity = 0.
     */
    StepCode
    (
     StepCode && other
     ) : core_(), bitCapacity_(other.bitCapacity_), bitSize_(other.bitSize_), size_(other.size_) {
      core_.mvStepCodeCore(std::move(other.core_), other.bitCapacity_, other.bitSize_, other.size_, bitCapacity_, bitSize_, size_);
      other.bitCapacity_ = other.bitSize_ = other.size_ = 0;
    }


    /*!
     * @brief Move operator.
     * @attention
     *   - 'other' is initialized to an object with capacity = 0.
     *   - The original contents are discarded.
     */
    StepCode & operator=
    (
     StepCode && other
     ) {
      if (this != &other) {
        core_.mvStepCodeCore(std::move(other.core_), other.bitCapacity_, other.bitSize_, other.size_, bitCapacity_, bitSize_, size_);
        other.bitCapacity_ = other.bitSize_ = other.size_ = 0;
      }
      return *this;
    }


    ~StepCode()
    {
    }


    /*!
     * @brief Get read-only array pointer.
     */
    const uint64_t * getConstPtr_vals() const noexcept
    {
      return core_.getConstPtr_vals();
    }


    /*!
     * @brief Get read-only wCodes_ array pointer.
     */
    const uint64_t * getConstPtr_wCodes() const noexcept
    {
      return core_.getConstPtr_wCodes();
    }


    uint64_t sumW
    (
     const uint64_t beg, //!< in [0, bv_.size()).
     const uint64_t end //!< in [beg, bv_.size()).
     ) const noexcept {
      assert(beg <= end);
      assert(end <= size_);

      return core_.sumW(beg, end);
    }


    /*!
     * @brief Read bit-width of "idx"-th (0base) values.
     */
    uint8_t readW
    (
     const uint64_t idx //!< in [0, size_).
     ) const noexcept {
      assert(idx < size_);

      return core_.readW(idx);
    }


    /*!
     * @brief Write wCode for "idx"-th (0base) value.
     */
    void writeWCode
    (
     const uint8_t wCode,
     const uint64_t idx //!< "(idx + 1) * 4 - 1"-th bit of "wCodes_" should not be out-of-bound.
     ) noexcept {
      assert(wCode <= ctcbits::UINTW_MAX(StepCodeUtil::kWCBits));

      core_.writeWCode(wCode, idx);
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in bv.
     */
    uint64_t calcBitPos
    (
     const uint64_t idx //!< "(idx + 1) * 4 - 1"-th bit of wCodesArray should not be out-of-bound.
     ) const noexcept {
      assert(idx < size_);

      return core_.calcBitPos(idx);
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in bv.
     */
    uint64_t calcBitPos
    (
     const uint64_t idx, //!< "(idx + 1) * 4 - 1"-th bit of wCodesArray should not be out-of-bound.
     const uint8_t * wCodesAuxM
     ) const noexcept {
      return core_.calcBitPos(idx, wCodesAuxM);
    }


    /*!
     * @brief Read "w"-bits written in the bit-region beginning at vals_[[bitPos..]] in bv_.
     * @return Value represented by vals_[[bitPos..bitPos+w)).
     * @pre The bit-region must not be out of bounds.
     */
    uint64_t readWBits
    (
     const uint64_t bitPos, //!< Bit-pos specifying the beginning position of the bit-region
     const uint8_t w //!< Bit-width in [0, 64].
     ) const noexcept {
      assert(bitPos < bitCapacity_);
      assert(bitPos + w <= bitCapacity_);
      assert(w <= 64);

      return core_.readWBits(bitPos, w);
    }


    /*!
     * @brief Write "w"-bit value "val" to the bit-region beginning at vals_[[bitPos..]].
     * @pre The bit-region must not be out of bounds.
     */
    void writeWBits
    (
     const uint64_t val, //!< in [0, 2^w).
     const uint64_t bitPos, //!< Bit-pos.
     const uint8_t steppedW //!< Bit-width in {4, 8, 12, 16, ..., 64}.
     ) noexcept {
      assert(bitPos < bitCapacity_);
      assert(bitPos + steppedW <= bitCapacity_);
      assert(steppedW <= 64);
      assert(steppedW % 4 == 0);
      assert(bits::bitSize(val) <= steppedW);

      core_.writeWBits(val, bitPos, steppedW);
    }


    /*!
     * @brief Read "idx"-th (0base) value.
     * @note It may be inefficient for sequencial access. Use readW and readWBits in that case.
     */
    uint64_t read
    (
     uint64_t idx //!< in [0, size_).
     ) const noexcept {
      assert(idx < size_);

      return core_.read(idx);
    }


    /*!
     * @brief Read "idx"-th (0base) value.
     * @note It may be inefficient for sequencial access. Use readW and readWBits in that case.
     */
    uint64_t read_naive
    (
     uint64_t idx //!< in [0, size_).
     ) const noexcept {
      assert(idx < size_);

      return core_.read_naive(idx);
    }


    /*!
     * @brief Write uint "val" at the end.
     * @pre The resulting size_ and bitSize_ should be within capacity.
     */
    void append
    (
     const uint64_t val, //!< Uint val.
     const uint8_t steppedW //!<  bit-width in {4, 8, 12, ..., 64} into which val fits.
     ) {
      assert(size_ < kMaxNum);
      assert(bitSize_ + steppedW <= bitCapacity_);
      assert(steppedW <= 64 && steppedW % 4 == 0);
      assert(bits::bitSize(val) <= steppedW);

      const uint8_t wCode = steppedW / StepCodeUtil::kStep - 1;
      core_.writeWBits(val, bitSize_, steppedW);
      core_.writeWCode(wCode, size_);
      bitSize_ += steppedW;
      ++size_;
    }


    /*!
     * @brief Write uint "val" at the end.
     * @pre The resulting size_ and bitSize_ should be within capacity.
     */
    void append
    (
     const uint64_t val //!< Uint val.
     ) {
      assert(size_ < kMaxNum);
      assert(bitSize_ + (bits::bitSize(val) + 3) / StepCodeUtil::kStep * StepCodeUtil::kStep <= bitCapacity_);

      this->append(val, StepCodeUtil::calcSteppedW(val));
    }



    /*!
     * @brief Rewrite value at "idx" to "val" without changing bit-width.
     */
    void rewriteVal
    (
     const uint64_t val, //!< should fit in readW(idx) bits.
     const uint64_t idx, //!< in [0, size_).
     const uint64_t bitPos //!< should be bit-position of "idx"-th (0base) values in bv_.
     ) {
      assert(idx < size_);
      assert(bits::bitSize(val) <= UINTW_MAX(readW(idx)));

      core_.rewriteVal(val, idx, bitPos);
    }


    /*!
     * @brief Get current capacity.
     * @pre The resulting size_ and bv_.size_ should be within capacity.
     */
    void changeWCodesAndValPos
    (
     const uint64_t * src,
     const uint64_t srcIdxBeg, //!< Beginning idx of src.
     const uint64_t srcLen, //!< Length of wCodes of src to insert.
     const uint64_t tgtIdxBeg, //!< Beginning idx of tgt.
     const uint64_t tgtLen, //!< Length of wCodes of tgt to delete.
     const uint64_t bitPos, //!< bitPos of bv_.
     const uint64_t insBitLen, //!< Bit-length to insert in bv_
     const uint64_t delBitLen //!< Bit-length to delete in bv_
     ) {
      assert(tgtIdxBeg + tgtLen <= size_);
      assert(size_ - tgtLen + srcLen < kMaxNum);
      assert(size_ + srcLen >= tgtLen);

      { // Change wCodes
        const auto tailNum = size_ - (tgtIdxBeg + tgtLen); // at least 0 by assumption.
        if (srcLen != tgtLen && tailNum) { // Need to shift elements in wCodes_.
          core_.mvWCodes(core_.getConstPtr_wCodes(), tgtIdxBeg + tgtLen, tgtIdxBeg + srcLen, tailNum);
        }
        if (srcLen) {
          core_.mvWCodes(src, srcIdxBeg, tgtIdxBeg, srcLen);
        }
        size_ += srcLen - tgtLen;
      }
      { // Shift positions of values
        const uint64_t srcPos = bitPos + delBitLen;
        const uint64_t tgtPos = bitPos + insBitLen;
        core_.mvVals(core_.getConstPtr_vals(), srcPos, tgtPos, bitSize_ - srcPos);
        bitSize_ += insBitLen - delBitLen;
      }
    }

  
    /*!
     * @brief Move vals.
     * @pre Bit-regions should be within capacity.
     */
    void mvVals
    (
     const uint64_t * srcVals,
     const uint64_t srcBitPos,
     const uint64_t tgtBitPos,
     const uint64_t bitLen
     ) noexcept {
      core_.mvVals(srcVals, srcBitPos, tgtBitPos, bitLen);
    }


    /*!
     * @brief Get current capacity of bv_.
     */
    size_t bitCapacity() const noexcept {
      return bitCapacity_;
    }


    /*!
     * @brief Get current size of bv_.
     */
    size_t bitSize() const noexcept {
      return bitSize_;
    }


    /*!
     * @brief Get current capacity.
     */
    size_t capacity() const noexcept {
      return kMaxNum;
    }


    /*!
     * @brief Get current size.
     */
    size_t size() const noexcept {
      return size_;
    }


    /*!
     * @brief Calculate total memory usage in bytes.
     */
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      return size + bitCapacity_ / 8;
    }


    /*!
     * @brief Return if no element is stored.
     */
    bool empty() const noexcept {
      return (size_ == 0);
    }


    /*!
     * @brief Clear vector. It only changes size to zero.
     */
    void clear() noexcept {
      size_ = 0;
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     * @node If givenCapacity is 0, it works as shrink_to_fit.
     */
    void changeBitCapacity
    (
     const size_t givenCapacity
     ) {
      assert(givenCapacity <= ctcbits::UINTW_MAX(sizeof(SizeT) * 8));
      assert(givenCapacity <= ctcbits::UINTW_MAX(58));

      if (bitCapacity_ != givenCapacity) {
        bitCapacity_ = static_cast<SizeT>(core_.setBitCapacity(std::max(bitSize_, static_cast<SizeT>(givenCapacity))));
      }
    }


    void printStatistics
    (
     const bool verbose = false
     ) const noexcept {
      std::cout << "StepCode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      std::cout << "size = " << this->size() << ", capacity = " << this->capacity() << std::endl;
      std::cout << this->calcMemBytes() << " bytes (dynamic array = " << sizeof(uint64_t) * (bitCapacity_ / 64) << ")" << std::endl;
      if (verbose) {
        {
          const auto size = this->size();
          std::cout << "dump bit witdth stored in wCodes (" << core_.getConstPtr_wCodes() << ")" << std::endl;
          for (uint64_t i = 0; i < size_; ++i) {
            std::cout << (uint64_t)(readW(i)) << " ";
          }
          std::cout << std::endl;
        }
        {
          const auto size = this->size();
          std::cout << "dump values" << std::endl;
          for (uint64_t i = 0; i < size_; ++i) {
            std::cout << read(i) << " ";
          }
          std::cout << std::endl;
        }
        {
          const auto size = this->bitSize();
          std::cout << "bitSize = " << size << ", bitCapacity = " << this->bitCapacity() << std::endl;
          std::cout << "dump bits in vals_ (" << getConstPtr_vals() << ")" << std::endl;
          for (uint64_t i = 0; i < (size + 63) / 64; ++i) {
            std::cout << "(" << i << ")";
            for (uint64_t j = 0; j < 64; ++j) {
              std::cout << bits::readWBits_S(getConstPtr_vals(), 64 * i + 63 - j, ctcbits::UINTW_MAX(1));
            }
            std::cout << " ";
          }
          std::cout << std::endl;
        }
      }
      std::cout << "StepCode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };
} // namespace itmmti

#endif

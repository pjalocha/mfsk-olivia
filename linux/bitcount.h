// Fast bit counting for (forward) error correction codes
// (c) 2003, Pawel Jalocha, Pawel.Jalocha@cern.ch

#ifndef __BITCOUNT_H__
#define __BITCOUNT_H__

#include <stdint.h>

// ==========================================================================
// a table for fast bit counting

const uint8_t ByteCount1s[256] = {
 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
} ;

// ==========================================================================

inline char Count1s(uint8_t Byte) { return ByteCount1s[Byte]; }

inline char Count1s(int8_t Byte) { return Count1s((uint8_t)Byte); }

inline char Count1s(uint16_t Word)
{ return Count1s((uint8_t)Word)
        +Count1s((uint8_t)(Word>>8)); }

inline char Count1s(int16_t Word) { return Count1s((uint16_t)Word); }

inline char Count1s(uint32_t LongWord)
{ return Count1s((uint16_t)LongWord)
        +Count1s((uint16_t)(LongWord>>16)); }

inline char Count1s(int32_t LongWord) { return Count1s((uint32_t)LongWord); }

// ==========================================================================

#endif // of __BITCOUNT_H__

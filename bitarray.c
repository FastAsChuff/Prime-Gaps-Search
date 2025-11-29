
//======================================================================================================================
// Bit Array: Functions for working with an array of booleans in c that only occupy 1 bit.
//======================================================================================================================
// Author: Simon Goater Nov 2025
// COPYRIGHT NOTICE: Copying, modifying, and distributing with conspicuous attribution for any purpose is permitted.
//======================================================================================================================
  #include<stdint.h>
  #include<stdbool.h>
  uint8_t *makebitarray(uint64_t numbits, uint64_t *numbytesret) {
    uint64_t numbytes = (numbits + 7U)/8U;
    uint8_t *retptr = calloc(1, numbytes);
    if (retptr) {
      if (numbytesret) *numbytesret = numbytes;
      return retptr;
    } else {
      return NULL;
    }
  }
  void setbitarray(uint64_t bitnum, uint8_t *bitarray) {
    uint64_t bytenum = bitnum / 8U;
    uint8_t bit = bitnum & 0x7ULL;
    uint8_t mask = 1U << (7U-bit);
    bitarray[bytenum] |= mask;
  }
  void unsetbitarray(uint64_t bitnum, uint8_t *bitarray) {
    uint64_t bytenum = bitnum / 8U;
    uint8_t bit = bitnum & 0x7ULL;
    uint8_t mask = 1U << (7U-bit);
    bitarray[bytenum] &= ~mask;
  }
  _Bool getbitarray(uint64_t bitnum, uint8_t *bitarray) {
    uint64_t bytenum = bitnum / 8U;
    uint8_t bit = bitnum & 0x7ULL;
    uint8_t mask = 1U << (7U-bit);
    return bitarray[bytenum] & mask;
  }
//======================================================================================================================


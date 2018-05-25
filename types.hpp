#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>
#include <cstddef>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t b8;
typedef uint16_t b16;
typedef uint32_t b32;
typedef uint64_t b64;

typedef float f32;
typedef double f64;
typedef long double f80;

typedef unsigned int uint;

#define ARR_SIZE(x) (sizeof(x)/sizeof(x[0]))

#include <cassert>

#define DEBUG_ASSERT assert

//define more descriptive keywords for the different uses of "static" keyword (thanks Casey!)
#define persist static
#define internal static

//NOTE: For now we are mostly only using "internal" for variable, because for functions in
//      header files, the "inline" keyword is probably better suited.
//NOTE: Some compilers (LLVM) do still interpret "inline" as an optimization hint,
//      which is not really something we want (we're using it only for linkage semantics).
//      However, we can replace the "inline" keyword using a global find-and-replace
//      without too many issues, so I'm not worried about changing my mind in the future.

#endif

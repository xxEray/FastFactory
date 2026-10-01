#ifndef _COMMON_
#define _COMMON_

#include <cstdio>
#include <cstdlib>

typedef unsigned long long ull;
typedef unsigned __int128 u128;

constexpr int MAX_BIT_FACTORY = 29;

constexpr int LEAF_PIVOT = 32;
constexpr int FAC_BF_PIVOT = 16;
constexpr int MUL_BF_PIVOT = 128;
constexpr int OMP_PIVOT = 4096;
// constexpr int SAVE_MEMORY = (1 << 27);
constexpr int LOG2_OMP_PIVOT = 12;
static_assert(1 << LOG2_OMP_PIVOT == OMP_PIVOT);

#ifdef DEBUG
#define assert(condition) ((condition) ? (void)0 : \
	(fprintf(stderr, "Assertion failed at ilne #%d in file " __FILE__ ": " #condition,\
	__LINE__), exit(-1), (void)0))
#else
#define assert(condition) ((void)0)
#endif

struct Trigger {
	template<typename T> Trigger(T &&func) { func(); }
};

#endif
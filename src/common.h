#ifndef _COMMON_
#define _COMMON_

#include <cstdio>
#include <cstdlib>

typedef unsigned long long ull;
typedef unsigned __int128 u128;

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
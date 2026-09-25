#ifndef _MULMOD_
#define _MULMOD_

#include "common.h"

const ull MOD1 = 1945555039024054273;
const ull G1 = 5;
const ull invG1 = 1167333023414432564;
static_assert((u128)G1 * invG1 % MOD1 == 1);

const ull MOD2 = 4179340454199820289;
const ull G2 = 3;
const ull invG2 = 1393113484733273430;
static_assert((u128)G2 * invG2 % MOD2 == 1);

const ull INV_MOD1_MOD2 = 3505253284167591212;
static_assert((u128)INV_MOD1_MOD2 * MOD1 % MOD2 == 1);

ull trim1(ull x) { return x >= MOD1 ? x - MOD1 : x;  }
void mulmod1(ull &x, ull y) {
	x = static_cast<u128>(x) * y % MOD1;
}

ull trim2(ull x) { return x >= MOD2 ? x - MOD2 : x; }
void mulmod2(ull &x, ull y) {
	x = static_cast<u128>(x) * y % MOD2;
}

#endif
#ifndef _MULMOD_
#define _MULMOD_

#include "common.h"

constexpr ull MOD1 = 1945555039024054273;
constexpr ull iMOD1 = 1945555039024054271; // MOD1 * iMOD1 = -1 (mod R = 2^64)
constexpr ull R_MOD1 = (static_cast<u128>(1) << 64) % MOD1;
constexpr ull R2_MOD1 = (static_cast<u128>(R_MOD1) << 64) % MOD1;
constexpr ull G1 = 5;
constexpr ull G1_R = (static_cast<u128>(G1) << 64) % MOD1;
constexpr ull invG1 = 1167333023414432564;
constexpr ull invG1_R = (static_cast<u128>(invG1) << 64) % MOD1;
static_assert(MOD1 * iMOD1 == static_cast<ull>(-1)); // (mod R = 2^64)
static_assert((u128)G1 * invG1 % MOD1 == 1);

constexpr ull MOD2 = 4179340454199820289;
constexpr ull iMOD2 = 4179340454199820287; // MOD2 * iMOD2 = -1 (mod R = 2^64)
constexpr ull R_MOD2 = (static_cast<u128>(1) << 64) % MOD2;
constexpr ull R2_MOD2 = (static_cast<u128>(R_MOD2) << 64) % MOD2;
constexpr ull G2 = 3;
constexpr ull G2_R = (static_cast<u128>(G2) << 64) % MOD2;
constexpr ull invG2 = 1393113484733273430;
constexpr ull invG2_R = (static_cast<u128>(invG2) << 64) % MOD2;
static_assert(MOD2 * iMOD2 == static_cast<ull>(-1)); // (mod R = 2^64)
static_assert((u128)G2 * invG2 % MOD2 == 1);

constexpr ull INV_MOD1_MOD2 = 3505253284167591212;
static_assert((u128)INV_MOD1_MOD2 * MOD1 % MOD2 == 1);

ull trim1(ull x) { return x >= MOD1 ? x - MOD1 : x;  }
void mulmod1(ull &x, ull y) {
	u128 v = static_cast<u128>(x) * y;
	ull t = static_cast<ull>(v) * iMOD1; // (mod R)
	x = (v + static_cast<u128>(MOD1) * t) >> 64;
	if(x >= MOD1) x -= MOD1;
}
void to_mont1(ull &x) { mulmod1(x, R2_MOD1); }
void from_mont1(ull &x) { mulmod1(x, 1); }

ull trim2(ull x) { return x >= MOD2 ? x - MOD2 : x; }
void mulmod2(ull &x, ull y) {
	u128 v = static_cast<u128>(x) * y;
	ull t = static_cast<ull>(v) * iMOD2; // (mod R)
	x = (v + static_cast<u128>(MOD2) * t) >> 64;
	if(x >= MOD2) x -= MOD2;
}
void to_mont2(ull &x) { mulmod2(x, R2_MOD2); }
void from_mont2(ull &x) { mulmod2(x, 1); }

#endif
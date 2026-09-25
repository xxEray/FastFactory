#ifndef _BIGINT_
#define _BIGINT_

#include <vector>
#include <algorithm>

#include "common.h"
#include "mulmod.h"

const ull B = 100000000000000ull;

ull qpow1(ull x, ull y) {
	ull ret = 1;
	for(; y; mulmod1(x, x), y >>= 1) if(y & 1) mulmod1(ret, x);
	return ret;
}
ull qpow2(ull x, ull y) {
	ull ret = 1;
	for(; y; mulmod2(x, x), y >>= 1) if(y & 1) mulmod2(ret, x);
	return ret;
}

struct BigInt {
	std::vector<ull> w;
	BigInt() : w({0}) {}
	BigInt(int val) : w({val}) {}
	BigInt(const BigInt &) = default;
	BigInt(BigInt &&) = default;
	BigInt &operator=(const BigInt &) = default;
	BigInt &operator=(BigInt &&) = default;
	void pop_zero() { while(w.size() > 1 && w.back() == 0) w.pop_back(); }
	BigInt &operator*=(const BigInt &rhs);
};

std::vector<int> rev[27];
void get_rev() {
	for(int w = 0; w < 27; w++) {
		rev[w].resize(1 << w);
		for(int i = 0; i < (1 << w); i++)
			rev[w][i] = w ? ((rev[w][i >> 1] >> 1) | ((i & 1) << (w - 1))) : 0;
	}
}
Trigger rev_trigger(get_rev);

BigInt operator*(const BigInt &x, const BigInt &y) {
	int len = x.w.size() + y.w.size();
	BigInt ans;
	if(len <= 8) {
		std::vector<u128> vec(len);
		for(int i = 0; i < (int)x.w.size(); i++)
			for(int j = 0; j < (int)y.w.size(); j++)
				vec[i + j] += static_cast<u128>(x.w[i]) * y.w[j];
		ans.w.resize(len);
		assert(len > 0);
		for(int i = 0; i < len - 1; i++) vec[i + 1] += vec[i] / B, ans.w[i] = vec[i] % B;
		ans.w[len - 1] = vec[len - 1] % B;
		ans.pop_zero();
	} else {
		int lim = 1, width = 0;
		while(lim < len) lim <<= 1, width++;
		std::vector<ull> vx = x.w, vy = y.w, ret1(lim), ret2(lim);
		vx.resize(lim), vy.resize(lim);
		auto NTT1 = [lim, width](std::vector<ull> &vec) {
			for(int i = 0; i < lim; i++)
				if(i < rev[width][i]) std::swap(vec[i], vec[rev[width][i]]);
			for(int i = 1; i < lim; i <<= 1) {
				ull wn = qpow1(G1, (MOD1 - 1) / (i << 1));
				for(int j = 0; j < lim; j += (i << 1)) {
					ull w = 1;
					for(int k = 0; k < i; k++, mulmod1(w, wn)) {
						ull x = vec[j + k], y = vec[j + i + k];
						mulmod1(y, w);
						vec[j + k] = trim1(x + y);
						vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
					}
				}
			}
		};
		NTT1(vx), NTT1(vy);
		for(int i = 0; i < lim; i++) ret1[i] = vx[i], mulmod1(ret1[i], vy[i]);
		auto iNTT1 = [lim, width](std::vector<ull> &vec) {
			for(int i = 0; i < lim; i++)
				if(i < rev[width][i]) std::swap(vec[i], vec[rev[width][i]]);
			for(int i = 1; i < lim; i <<= 1) {
				ull wn = qpow1(invG1, (MOD1 - 1) / (i << 1));
				for(int j = 0; j < lim; j += (i << 1)) {
					ull w = 1;
					for(int k = 0; k < i; k++, mulmod1(w, wn)) {
						ull x = vec[j + k], y = vec[j + i + k];
						mulmod1(y, w);
						vec[j + k] = trim1(x + y);
						vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
					}
				}
			}
			ull inv = qpow1(lim, MOD1 - 2);
			for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
		};
		iNTT1(ret1);
		auto NTT2 = [lim, width](std::vector<ull> &vec) {
			for(int i = 0; i < lim; i++)
				if(i < rev[width][i]) std::swap(vec[i], vec[rev[width][i]]);
			for(int i = 1; i < lim; i <<= 1) {
				ull wn = qpow2(G2, (MOD2 - 1) / (i << 1));
				for(int j = 0; j < lim; j += (i << 1)) {
					ull w = 1;
					for(int k = 0; k < i; k++, mulmod2(w, wn)) {
						ull x = vec[j + k], y = vec[j + i + k];
						mulmod2(y, w);
						vec[j + k] = trim2(x + y);
						vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
					}
				}
			}
		};
		vx = x.w, vy = y.w;
		vx.resize(lim), vy.resize(lim);
		NTT2(vx), NTT2(vy);
		for(int i = 0; i < lim; i++) ret2[i] = vx[i], mulmod2(ret2[i], vy[i]);
		auto iNTT2 = [lim, width](std::vector<ull> &vec) {
			for(int i = 0; i < lim; i++)
				if(i < rev[width][i]) std::swap(vec[i], vec[rev[width][i]]);
			for(int i = 1; i < lim; i <<= 1) {
				ull wn = qpow2(invG2, (MOD2 - 1) / (i << 1));
				for(int j = 0; j < lim; j += (i << 1)) {
					ull w = 1;
					for(int k = 0; k < i; k++, mulmod2(w, wn)) {
						ull x = vec[j + k], y = vec[j + i + k];
						mulmod2(y, w);
						vec[j + k] = trim2(x + y);
						vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
					}
				}
			}
			ull inv = qpow2(lim, MOD2 - 2);
			for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
		};
		iNTT2(ret2);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		ull last = 0;
		for(int i = 0; i < lim; i++) {
			ull t = (ret2[i] >= ret1[i] ? ret2[i] - ret1[i] : ret2[i] + MOD2 - ret1[i]);
			mulmod2(t, INV_MOD1_MOD2);
			u128 val = static_cast<u128>(ret1[i]) + static_cast<u128>(MOD1) * t;
			val += last;
			last = val / B;
			ret1[i] = static_cast<ull>(val % B);
		}
		ans.w = ret1;
		ans.w.resize(len);
		ans.pop_zero();
	}
	return ans;
}

inline BigInt &BigInt::operator*=(const BigInt &rhs) { return *this = *this * rhs; }

#endif
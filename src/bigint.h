#ifndef _BIGINT_
#define _BIGINT_

#include <vector>
#include <algorithm>

#include "common.h"
#include "speed.h"
#include "mulmod.h"

constexpr ull B = 100000000000000ull; // 1e14

ull qpow1(ull x, ull y) {
	ull ret = R_MOD1; // to_mont1(1)
	for(; y; mulmod1(x, x), y >>= 1) if(y & 1) mulmod1(ret, x);
	return ret;
}
ull qpow2(ull x, ull y) {
	ull ret = R_MOD2; // to_mont2(1)
	for(; y; mulmod2(x, x), y >>= 1) if(y & 1) mulmod2(ret, x);
	return ret;
}

struct BigInt {
	std::vector<ull> w;
	BigInt() : w({0}) {}
	BigInt(unsigned val) : w({val}) {}
	BigInt(const BigInt &) = default;
	BigInt(BigInt &&) = default;
	BigInt &operator=(const BigInt &) = default;
	BigInt &operator=(BigInt &&) = default;
	void pop_zero() { while(w.size() > 1 && w.back() == 0) w.pop_back(); }
	BigInt &operator*=(unsigned rhs) {
		w.emplace_back(0);
		ull last = 0;
		for(int i = 0; i < (int)w.size(); i++) {
			u128 val = static_cast<u128>(w[i]) * rhs + last;
			w[i] = val % B, last = val / B;
		}
		pop_zero();
		return *this;
	}
};

ull invlim1[MAX_BIT_FACTORY], invlim2[MAX_BIT_FACTORY];
ull wn1[MAX_BIT_FACTORY], iwn1[MAX_BIT_FACTORY], wn2[MAX_BIT_FACTORY], iwn2[MAX_BIT_FACTORY];
ull tw1[LOG2_OMP_PIVOT][OMP_PIVOT], tw2[LOG2_OMP_PIVOT][OMP_PIVOT];
std::vector<ull> tmp1, tmp2;

void get_wn() {
	for(int w = 0; w < MAX_BIT_FACTORY; w++) {
		// rev[w].resize(1 << w);
		// for(int i = 0; i < (1 << w); i++)
		// 	rev[w][i] = w ? ((rev[w][i >> 1] >> 1) | ((i & 1) << (w - 1))) : 0;
		ull lim1 = (1 << w), lim2 = (1 << w);
		to_mont1(lim1), to_mont2(lim2);
		invlim1[w] = qpow1(lim1, MOD1 - 2), invlim2[w] = qpow2(lim2, MOD2 - 2);
		wn1[w] = qpow1(G1_R, (MOD1 - 1) / (1 << w));
		iwn1[w] = qpow1(invG1_R, (MOD1 - 1) / (1 << w));
		wn2[w] = qpow2(G2_R, (MOD2 - 1) / (1 << w));
		iwn2[w] = qpow2(invG2_R, (MOD2 - 1) / (1 << w));
	}
	for(int w = 0; w < LOG2_OMP_PIVOT; w++) {
		auto calc1 = [w](ull vec[OMP_PIVOT], ull wn) {
			vec[0] = R_MOD1;
			for(int i = 1; i < OMP_PIVOT; i++) vec[i] = vec[i - 1], mulmod1(vec[i], wn);
		};
		auto calc2 = [w](ull vec[OMP_PIVOT], ull wn) {
			vec[0] = R_MOD2;
			for(int i = 1; i < OMP_PIVOT; i++) vec[i] = vec[i - 1], mulmod2(vec[i], wn);
		};
		calc1(tw1[w], wn1[w]);
		calc2(tw2[w], wn2[w]);
	}
}
Trigger rev_trigger(get_wn);

void mul_eq(BigInt &x, BigInt &&y) {
	int len = x.w.size() + y.w.size();
	if(len <= MUL_BF_PIVOT) {
#ifdef SPEED
		++speed_rec[0].calls;
#endif
		BigInt ans;
		std::vector<u128> vec(len);
		for(int i = 0; i < (int)x.w.size(); i++)
			for(int j = 0; j < (int)y.w.size(); j++)
				vec[i + j] += static_cast<u128>(x.w[i]) * y.w[j];
		ans.w.resize(len);
		assert(len > 0);
		for(int i = 0; i < len - 1; i++) vec[i + 1] += vec[i] / B, ans.w[i] = vec[i] % B;
		ans.w[len - 1] = vec[len - 1] % B;
		ans.pop_zero();
		std::swap(ans.w, x.w);
	} else {
		int lim = 1, width = 0;
		while(lim < len) lim <<= 1, width++;
		assert(width < MAX_BIT_FACTORY);
		x.w.resize(lim), y.w.resize(lim), tmp1 = x.w;
		/* if(lim < SAVE_MEMORY)  */tmp2 = y.w;
		auto &v1 = x.w, &v2 = y.w, &v3 = tmp1, &v4 = tmp2;
#ifdef SPEED
		SpeedRec &sr = speed_rec[lim];
		++sr.calls;
#endif
		if(lim < OMP_PIVOT) {
			auto DIF1 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						// ull wn = wn1[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							// ull w = R_MOD1; // to_mont1(1)
							for(int k = 0; k < i; k++/* , mulmod1(w, wn) */) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
								mulmod1(vec[j + i + k], tw1[cn][k]);
							}
						}
					}
				}
			};
			SPEED_TICK(a0);
			for(int i = 0; i < lim; i++) to_mont1(v3[i]), to_mont1(v4[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			DIF1(v3), DIF1(v4);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) mulmod1(v3[i], v4[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			auto DIT1 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn1[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1;
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod1(y, w);
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
							}
						}
					}
				}
				ull inv = invlim1[width];
				for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
			};
			SPEED_TICK(a5);
			DIT1(v3);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);
			auto DIF2 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						// ull wn = wn2[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							// ull w = R_MOD2;
							for(int k = 0; k < i; k++/* , mulmod2(w, wn) */) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
								mulmod2(vec[j + i + k], tw2[cn][k]);
							}
						}
					}
				}
			};
			SPEED_TICK(b0);
			for(int i = 0; i < lim; i++) to_mont2(v1[i]), to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			DIF2(v1), DIF2(v2);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) mulmod2(v1[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			auto DIT2 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn2[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod2(y, w);
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
							}
						}
					}
				}
				ull inv = invlim2[width];
				for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
			};
			SPEED_TICK(b5);
			DIT2(v1);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		} else /* if(lim < SAVE_MEMORY)  */{
			auto DIF1 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1; // to_mont1(1)
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
								mulmod1(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(a0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) to_mont1(v3[i]), to_mont1(v4[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			DIF1(v3), DIF1(v4);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod1(v3[i], v4[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			auto DIT1 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1;
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod1(y, w);
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
							}
						}
					}
				}
				ull inv = invlim1[width];
				#pragma omp parallel for
				for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
			};
			SPEED_TICK(a5);
			DIT1(v3);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);
			auto DIF2 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
								mulmod2(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(b0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) to_mont2(v1[i]), to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			DIF2(v1), DIF2(v2);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod2(v1[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			auto DIT2 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod2(y, w);
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
							}
						}
					}
				}
				ull inv = invlim2[width];
				#pragma omp parallel for
				for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
			};
			SPEED_TICK(b5);
			DIT2(v1);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}/*  else {
			auto DIF1 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1; // to_mont1(1)
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
								mulmod1(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(a0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) to_mont1(v3[i]), to_mont1(v2[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			DIF1(v3), DIF1(v2);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod1(v3[i], v2[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			auto DIT1 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1;
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod1(y, w);
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
							}
						}
					}
				}
				ull inv = invlim1[width];
				#pragma omp parallel for
				for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
			};
			SPEED_TICK(a5);
			DIT1(v3), DIT1(v2);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);
			auto DIF2 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
								mulmod2(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(b0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++)
				from_mont1(v2[i]), to_mont2(v1[i]), to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			DIF2(v1), DIF2(v2);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod2(v1[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			auto DIT2 = [lim, width](std::vector<ull> &vec) {
				#pragma omp parallel
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod2(y, w);
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
							}
						}
					}
				}
				ull inv = invlim2[width];
				#pragma omp parallel for
				for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
			};
			SPEED_TICK(b5);
			DIT2(v1);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		} */
		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		u128 last = 0;
		for(int i = 0; i < lim; i++) {
			from_mont1(v3[i]), from_mont2(v1[i]);
			ull t = (v1[i] >= v3[i] ? v1[i] - v3[i] : v1[i] + MOD2 - v3[i]);
			// t = static_cast<u128>(t) * INV_MOD1_MOD2 % MOD2;
			mulmod2(t, INV_MOD1_MOD2_R);
			u128 val = static_cast<u128>(v3[i]) + static_cast<u128>(MOD1) * t;
			val += last;
			last = val / B;
			v1[i] = static_cast<ull>(val % B);
		}
		SPEED_TICK(t9);
		SPEED_ADD(crt, t8, t9);
		x.pop_zero();
	}
}

void mul_self_eq(BigInt &x) {
	int len = x.w.size() * 2;
	if(len <= MUL_BF_PIVOT) {
#ifdef SPEED
		++speed_rec_self[0].calls;
#endif
		BigInt ans;
		std::vector<u128> vec(len);
		for(int i = 0; i < (int)x.w.size(); i++)
			for(int j = 0; j < (int)x.w.size(); j++)
				vec[i + j] += static_cast<u128>(x.w[i]) * x.w[j];
		ans.w.resize(len);
		assert(len > 0);
		for(int i = 0; i < len - 1; i++) vec[i + 1] += vec[i] / B, ans.w[i] = vec[i] % B;
		ans.w[len - 1] = vec[len - 1] % B;
		ans.pop_zero();
		std::swap(ans.w, x.w);
	} else {
		int lim = 1, width = 0;
		while(lim < len) lim <<= 1, width++;
		x.w.resize(lim), tmp1 = x.w;
		std::vector<ull> &v1 = x.w, &v2 = tmp1;
#ifdef SPEED
		SpeedRec &sr = speed_rec_self[lim];
		++sr.calls;
#endif
		if(lim < OMP_PIVOT) {
			auto DIF1 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						// ull wn = wn1[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							// ull w = R_MOD1; // to_mont1(1)
							for(int k = 0; k < i; k++/* , mulmod1(w, wn) */) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
								mulmod1(vec[j + i + k], tw1[cn][k]);
							}
						}
					}
				}
			};
			SPEED_TICK(a0);
			for(int i = 0; i < lim; i++) to_mont1(v1[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			DIF1(v1);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) mulmod1(v1[i], v1[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			auto DIT1 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn1[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1;
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod1(y, w);
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
							}
						}
					}
				}
				ull inv = invlim1[width];
				for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
			};
			SPEED_TICK(a5);
			DIT1(v1);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);
			auto DIF2 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						// ull wn = wn2[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							// ull w = R_MOD2;
							for(int k = 0; k < i; k++/* , mulmod2(w, wn) */) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
								mulmod2(vec[j + i + k], tw2[cn][k]);
							}
						}
					}
				}
			};
			SPEED_TICK(b0);
			for(int i = 0; i < lim; i++) to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			DIF2(v2);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) mulmod2(v2[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			auto DIT2 = [lim, width](std::vector<ull> &vec) {
				{
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn2[cn];
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod2(y, w);
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
							}
						}
					}
				}
				ull inv = invlim2[width];
				for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
			};
			SPEED_TICK(b5);
			DIT2(v2);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		} else {
			auto DIF1 = [lim, width](std::vector<ull> &vec) {
				{
					#pragma omp parallel
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1; // to_mont1(1)
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
								mulmod1(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(a0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) to_mont1(v1[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			DIF1(v1);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod1(v1[i], v1[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			auto DIT1 = [lim, width](std::vector<ull> &vec) {
				{
					#pragma omp parallel
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn1[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD1;
							for(int k = 0; k < i; k++, mulmod1(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod1(y, w);
								vec[j + k] = trim1(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD1 - y);
							}
						}
					}
				}
				ull inv = invlim1[width];
				for(int i = 0; i < lim; i++) mulmod1(vec[i], inv);
			};
			SPEED_TICK(a5);
			DIT1(v1);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);
			auto DIF2 = [lim, width](std::vector<ull> &vec) {
				{
					#pragma omp parallel
					for(int i = (lim >> 1), cn = width; i >= 1; i >>= 1, cn--) {
						ull wn = wn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
								mulmod2(vec[j + i + k], w);
							}
						}
					}
				}
			};
			SPEED_TICK(b0);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			DIF2(v2);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			#pragma omp parallel for
			for(int i = 0; i < lim; i++) mulmod2(v2[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			auto DIT2 = [lim, width](std::vector<ull> &vec) {
				{
					#pragma omp parallel
					for(int i = 1, cn = 1; i < lim; i <<= 1, cn++) {
						ull wn = iwn2[cn];
						#pragma omp for
						for(int j = 0; j < lim; j += (i << 1)) {
							ull w = R_MOD2;
							for(int k = 0; k < i; k++, mulmod2(w, wn)) {
								ull x = vec[j + k], y = vec[j + i + k];
								mulmod2(y, w);
								vec[j + k] = trim2(x + y);
								vec[j + i + k] = (x >= y ? x - y : x + MOD2 - y);
							}
						}
					}
				}
				ull inv = invlim2[width];
				for(int i = 0; i < lim; i++) mulmod2(vec[i], inv);
			};
			SPEED_TICK(b5);
			DIT2(v2);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}
		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		u128 last = 0;
		for(int i = 0; i < lim; i++) {
			from_mont1(v1[i]), from_mont2(v2[i]);
			ull t = (v2[i] >= v1[i] ? v2[i] - v1[i] : v2[i] + MOD2 - v1[i]);
			// t = static_cast<u128>(t) * INV_MOD1_MOD2 % MOD2;
			mulmod2(t, INV_MOD1_MOD2_R);
			u128 val = static_cast<u128>(v1[i]) + static_cast<u128>(MOD1) * t;
			val += last;
			last = val / B;
			v1[i] = static_cast<ull>(val % B);
		}
		SPEED_TICK(t9);
		SPEED_ADD(crt, t8, t9);
		x.pop_zero();
	}
}

#endif
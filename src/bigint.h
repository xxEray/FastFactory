#ifndef _BIGINT_
#define _BIGINT_

#include <vector>
#include <algorithm>

#include "common.h"
#include "speed.h"
#include "mulmod.h"

#ifdef SPEED
#include <chrono>
// Records the elapsed time since the previous mark into `field`; evaluated by
// the master thread only, so it is safe inside an OpenMP parallel region.
#define SPEED_MARK(field) do { \
	auto _t = std::chrono::steady_clock::now(); \
	SPEED_ADD(field, tprev, _t); \
	tprev = _t; \
} while(0)
#else
#define SPEED_MARK(field) ((void)0)
#endif

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
ull invlim1f[MAX_BIT_FACTORY], invlim2f[MAX_BIT_FACTORY]; // invlim with the Montgomery factor removed
ull wn1[MAX_BIT_FACTORY], iwn1[MAX_BIT_FACTORY], wn2[MAX_BIT_FACTORY], iwn2[MAX_BIT_FACTORY];
ull im1, im2;        // sqrt(-1) in Montgomery form
ull im1inv, im2inv;  // -sqrt(-1) in Montgomery form
// Pow tables (radix-2 / radix-4 twiddles) for short stage twiddle sequences.
ull tw1[LOG2_OMP_PIVOT][OMP_PIVOT], tw2[LOG2_OMP_PIVOT][OMP_PIVOT];
ull itw1[LOG2_OMP_PIVOT][OMP_PIVOT], itw2[LOG2_OMP_PIVOT][OMP_PIVOT];
ull twim1[LOG2_OMP_PIVOT][OMP_PIVOT], twim2[LOG2_OMP_PIVOT][OMP_PIVOT];
ull itwim1[LOG2_OMP_PIVOT][OMP_PIVOT], itwim2[LOG2_OMP_PIVOT][OMP_PIVOT];
// Scratch buffers. thread_local so that several independent multiplications
// may run concurrently (e.g. independent product-tree nodes on different
// threads) without clobbering each other.
thread_local std::vector<ull> tmp1, tmp2;

void get_wn() {
	for(int w = 0; w < MAX_BIT_FACTORY; w++) {
		ull lim1 = (1 << w), lim2 = (1 << w);
		to_mont1(lim1), to_mont2(lim2);
		invlim1[w] = qpow1(lim1, MOD1 - 2), invlim2[w] = qpow2(lim2, MOD2 - 2);
		invlim1f[w] = invlim1[w], from_mont1(invlim1f[w]);
		invlim2f[w] = invlim2[w], from_mont2(invlim2f[w]);
		wn1[w] = qpow1(G1_R, (MOD1 - 1) / (1 << w));
		iwn1[w] = qpow1(invG1_R, (MOD1 - 1) / (1 << w));
		wn2[w] = qpow2(G2_R, (MOD2 - 1) / (1 << w));
		iwn2[w] = qpow2(invG2_R, (MOD2 - 1) / (1 << w));
	}
	// sqrt(-1) = g^((p-1)/4); both moduli are 1 (mod 4).
	im1 = qpow1(G1_R, (MOD1 - 1) / 4), im1inv = MOD1 - im1;
	im2 = qpow2(G2_R, (MOD2 - 1) / 4), im2inv = MOD2 - im2;
	for(int w = 0; w < LOG2_OMP_PIVOT; w++) {
		tw1[w][0] = R_MOD1, itw1[w][0] = R_MOD1;
		tw2[w][0] = R_MOD2, itw2[w][0] = R_MOD2;
		for(int i = 1; i < OMP_PIVOT; i++) {
			tw1[w][i] = tw1[w][i - 1], mulmod1(tw1[w][i], wn1[w]);
			itw1[w][i] = itw1[w][i - 1], mulmod1(itw1[w][i], iwn1[w]);
			tw2[w][i] = tw2[w][i - 1], mulmod2(tw2[w][i], wn2[w]);
			itw2[w][i] = itw2[w][i - 1], mulmod2(itw2[w][i], iwn2[w]);
		}
		for(int i = 0; i < OMP_PIVOT; i++) {
			twim1[w][i] = tw1[w][i], mulmod1(twim1[w][i], im1);
			itwim1[w][i] = itw1[w][i], mulmod1(itwim1[w][i], im1inv);
			twim2[w][i] = tw2[w][i], mulmod2(twim2[w][i], im2);
			itwim2[w][i] = itw2[w][i], mulmod2(itwim2[w][i], im2inv);
		}
	}
}
Trigger rev_trigger(get_wn);

// ---- modulus traits -------------------------------------------------------
template<int ID> struct Mod;
template<> struct Mod<1> {
	static constexpr ull MOD = MOD1;
	static constexpr ull R = R_MOD1;
	static inline ull trim(ull x) { return x >= MOD1 ? x - MOD1 : x; }
	static inline void mulmod(ull &x, ull y) { mulmod1(x, y); }
	static inline ull wn(int c) { return wn1[c]; }
	static inline ull iwn(int c) { return iwn1[c]; }
	static inline ull im() { return im1; }
	static inline ull iminv() { return im1inv; }
	static inline ull *tw(int c) { return tw1[c]; }
	static inline ull *itw(int c) { return itw1[c]; }
	static inline ull *twim(int c) { return twim1[c]; }
	static inline ull *itwim(int c) { return itwim1[c]; }
	static inline ull invlimf(int w) { return invlim1f[w]; }
};
template<> struct Mod<2> {
	static constexpr ull MOD = MOD2;
	static constexpr ull R = R_MOD2;
	static inline ull trim(ull x) { return x >= MOD2 ? x - MOD2 : x; }
	static inline void mulmod(ull &x, ull y) { mulmod2(x, y); }
	static inline ull wn(int c) { return wn2[c]; }
	static inline ull iwn(int c) { return iwn2[c]; }
	static inline ull im() { return im2; }
	static inline ull iminv() { return im2inv; }
	static inline ull *tw(int c) { return tw2[c]; }
	static inline ull *itw(int c) { return itw2[c]; }
	static inline ull *twim(int c) { return twim2[c]; }
	static inline ull *itwim(int c) { return itwim2[c]; }
	static inline ull invlimf(int w) { return invlim2f[w]; }
};
using M1 = Mod<1>;
using M2 = Mod<2>;

// ---- radix-2/radix-4 butterflies ------------------------------------------
// A stage with parameter cn has block size 2^cn (half i = 2^(cn-1)).
// A radix-4 step with parameter cn fuses stages cn and cn-1 (block 4m, m = 2^(cn-2)).
// DIF (forward, decimation in frequency), no bit reversal.

// Shared radix-4 butterflies. Twiddles are passed in (Montgomery form):
//   DIF: p = W^k,     q = W^(2k),   r = W^k * IM
//   DIT: p = W^(-2k), q = W^(-k),   r = W^(-k) * IM^(-1)
template<class M>
static inline void dif4_butterfly(std::vector<ull> &vec, int o, int m, ull p, ull q, ull r) {
	ull a0 = vec[o], a1 = vec[o + m], a2 = vec[o + 2 * m], a3 = vec[o + 3 * m];
	ull t0 = M::trim(a0 + a2), t1 = (a0 >= a2 ? a0 - a2 : a0 + M::MOD - a2);
	ull t2 = M::trim(a1 + a3), t3 = (a1 >= a3 ? a1 - a3 : a1 + M::MOD - a3);
	ull c0 = M::trim(t0 + t2), c1 = (t0 >= t2 ? t0 - t2 : t0 + M::MOD - t2);
	ull A = t1, B = t3;
	M::mulmod(A, p);
	M::mulmod(B, r);
	ull c2 = M::trim(A + B), c3 = (A >= B ? A - B : A + M::MOD - B);
	M::mulmod(c1, q);
	M::mulmod(c3, q);
	vec[o] = c0, vec[o + m] = c1, vec[o + 2 * m] = c2, vec[o + 3 * m] = c3;
}
template<class M>
static inline void dit4_butterfly(std::vector<ull> &vec, int o, int m, ull p, ull q, ull r) {
	ull c0 = vec[o], c1 = vec[o + m], c2 = vec[o + 2 * m], c3 = vec[o + 3 * m];
	ull F = c1, E = c3;
	M::mulmod(F, p);
	M::mulmod(E, p);
	ull t0 = M::trim(c0 + F), t2 = (c0 >= F ? c0 - F : c0 + M::MOD - F);
	ull s = M::trim(c2 + E), d = (c2 >= E ? c2 - E : c2 + M::MOD - E);
	ull t1 = s, t3 = d;
	M::mulmod(t1, q);
	M::mulmod(t3, r);
	ull a0 = M::trim(t0 + t1), a2 = (t0 >= t1 ? t0 - t1 : t0 + M::MOD - t1);
	ull a1 = M::trim(t2 + t3), a3 = (t2 >= t3 ? t2 - t3 : t2 + M::MOD - t3);
	vec[o] = a0, vec[o + m] = a1, vec[o + 2 * m] = a2, vec[o + 3 * m] = a3;
}

template<class M>
void dif2_ser(std::vector<ull> &vec, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	const ull *twc = M::tw(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < i; k++) {
			ull x = vec[j + k], y = vec[j + i + k];
			vec[j + k] = M::trim(x + y);
			vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
			M::mulmod(vec[j + i + k], twc[k]);
		}
}

template<class M>
void dif4_ser(std::vector<ull> &vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	const ull *twc = M::tw(cn), *twp = M::tw(cn - 1), *twic = M::twim(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < m; k++)
			dif4_butterfly<M>(vec, j + k, m, twc[k], twp[k], twic[k]);
}

template<class M>
void dit2_ser(std::vector<ull> &vec, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	const ull *itwc = M::itw(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < i; k++) {
			ull x = vec[j + k], y = vec[j + i + k];
			M::mulmod(y, itwc[k]);
			vec[j + k] = M::trim(x + y);
			vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
		}
}

template<class M>
void dit4_ser(std::vector<ull> &vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	const ull *itwc = M::itw(cn), *itwp = M::itw(cn - 1), *itwic = M::itwim(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < m; k++)
			dit4_butterfly<M>(vec, j + k, m, itwp[k], itwc[k], itwic[k]);
}

// ---- OpenMP stages --------------------------------------------------------
template<class M>
void dif2_omp(std::vector<ull> &vec, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < LOG2_OMP_PIVOT) {
		const ull *twc = M::tw(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step)
			for(int k = 0; k < i; k++) {
				ull x = vec[j + k], y = vec[j + i + k];
				vec[j + k] = M::trim(x + y);
				vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
				M::mulmod(vec[j + i + k], twc[k]);
			}
	} else {
		ull wn = M::wn(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step) {
			ull w = M::R;
			for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
				ull x = vec[j + k], y = vec[j + i + k];
				vec[j + k] = M::trim(x + y);
				vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
				M::mulmod(vec[j + i + k], w);
			}
		}
	}
}

template<class M>
void dif4_omp(std::vector<ull> &vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < LOG2_OMP_PIVOT) {
		const ull *twc = M::tw(cn), *twp = M::tw(cn - 1), *twic = M::twim(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step)
			for(int k = 0; k < m; k++)
				dif4_butterfly<M>(vec, j + k, m, twc[k], twp[k], twic[k]);
	} else {
		ull W = M::wn(cn), W2 = M::wn(cn - 1), IM = M::im();
		#pragma omp for
		for(int j = 0; j < lim; j += step) {
			ull p = M::R, q = M::R, r = IM;
			for(int k = 0; k < m; k++) {
				dif4_butterfly<M>(vec, j + k, m, p, q, r);
				M::mulmod(p, W), M::mulmod(q, W2), M::mulmod(r, W);
			}
		}
	}
}

template<class M>
void dit2_omp(std::vector<ull> &vec, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < LOG2_OMP_PIVOT) {
		const ull *itwc = M::itw(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step)
			for(int k = 0; k < i; k++) {
				ull x = vec[j + k], y = vec[j + i + k];
				M::mulmod(y, itwc[k]);
				vec[j + k] = M::trim(x + y);
				vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
			}
	} else {
		ull wn = M::iwn(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step) {
			ull w = M::R;
			for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
				ull x = vec[j + k], y = vec[j + i + k];
				M::mulmod(y, w);
				vec[j + k] = M::trim(x + y);
				vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
			}
		}
	}
}

template<class M>
void dit4_omp(std::vector<ull> &vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < LOG2_OMP_PIVOT) {
		const ull *itwc = M::itw(cn), *itwp = M::itw(cn - 1), *itwic = M::itwim(cn);
		#pragma omp for
		for(int j = 0; j < lim; j += step)
			for(int k = 0; k < m; k++)
				dit4_butterfly<M>(vec, j + k, m, itwp[k], itwc[k], itwic[k]);
	} else {
		ull iW = M::iwn(cn), iW2 = M::iwn(cn - 1), IMI = M::iminv();
		#pragma omp for
		for(int j = 0; j < lim; j += step) {
			ull p = M::R, q = M::R, r = IMI; // p = W^(-2k), q = W^(-k), r = W^(-k) * IM^(-1)
			for(int k = 0; k < m; k++) {
				dit4_butterfly<M>(vec, j + k, m, p, q, r);
				M::mulmod(p, iW2), M::mulmod(q, iW), M::mulmod(r, iW);
			}
		}
	}
}

// ---- transform drivers ----------------------------------------------------
template<class M>
void dif_ser(std::vector<ull> &vec, int lim, int width) {
	int cn = width;
	if(width & 1) { dif2_ser<M>(vec, lim, cn); cn--; }
	for(; cn >= 2; cn -= 2) dif4_ser<M>(vec, lim, cn);
}
template<class M>
void dit_ser(std::vector<ull> &vec, int lim, int width) {
	int cmax = (width & 1) ? width - 1 : width;
	for(int cn = 2; cn <= cmax; cn += 2) dit4_ser<M>(vec, lim, cn);
	if(width & 1) dit2_ser<M>(vec, lim, width);
	ull inv = M::invlimf(width);
	for(int i = 0; i < lim; i++) M::mulmod(vec[i], inv);
}
// The *_stages drivers below contain `#pragma omp for` worksharing loops and
// MUST be called from inside an already-active `#pragma omp parallel` region
// (this lets one parallel region cover a whole multiplication instead of one
// region per array/transform).
template<class M>
void dif_stages(std::vector<ull> &vec, int lim, int width) {
	int cn = width;
	if(width & 1) { dif2_omp<M>(vec, lim, cn); cn--; }
	for(; cn >= 2; cn -= 2) dif4_omp<M>(vec, lim, cn);
}
template<class M>
void dit_stages(std::vector<ull> &vec, int lim, int width) {
	int cmax = (width & 1) ? width - 1 : width;
	for(int cn = 2; cn <= cmax; cn += 2) dit4_omp<M>(vec, lim, cn);
	if(width & 1) dit2_omp<M>(vec, lim, width);
	ull inv = M::invlimf(width);
	#pragma omp for
	for(int i = 0; i < lim; i++) M::mulmod(vec[i], inv);
}

// ---- multi-array transforms ------------------------------------------------
// Same stages as above, but a single worksharing loop covers several arrays at
// once (with `collapse(2)` over {array, block}).  This widens the parallelism
// of the early DIF / late DIT stages, which have only a handful of blocks each
// and would otherwise leave most threads idle.
struct BufRef { std::vector<ull> *v; bool m2; };

template<class M>
static inline void dif2_block(std::vector<ull> &vec, int j, int i, const ull *twc) {
	for(int k = 0; k < i; k++) {
		ull x = vec[j + k], y = vec[j + i + k];
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
		M::mulmod(vec[j + i + k], twc[k]);
	}
}
template<class M>
static inline void dif2_block_run(std::vector<ull> &vec, int j, int i, ull wn) {
	ull w = M::R;
	for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
		ull x = vec[j + k], y = vec[j + i + k];
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
		M::mulmod(vec[j + i + k], w);
	}
}
template<class M>
static inline void dif4_block(std::vector<ull> &vec, int j, int m, const ull *twc, const ull *twp, const ull *twic) {
	for(int k = 0; k < m; k++)
		dif4_butterfly<M>(vec, j + k, m, twc[k], twp[k], twic[k]);
}
template<class M>
static inline void dif4_block_run(std::vector<ull> &vec, int j, int m, ull W, ull W2, ull IM) {
	ull p = M::R, q = M::R, r = IM;
	for(int k = 0; k < m; k++) {
		dif4_butterfly<M>(vec, j + k, m, p, q, r);
		M::mulmod(p, W), M::mulmod(q, W2), M::mulmod(r, W);
	}
}
template<class M>
static inline void dit2_block(std::vector<ull> &vec, int j, int i, const ull *itwc) {
	for(int k = 0; k < i; k++) {
		ull x = vec[j + k], y = vec[j + i + k];
		M::mulmod(y, itwc[k]);
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
	}
}
template<class M>
static inline void dit2_block_run(std::vector<ull> &vec, int j, int i, ull wn) {
	ull w = M::R;
	for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
		ull x = vec[j + k], y = vec[j + i + k];
		M::mulmod(y, w);
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
	}
}
template<class M>
static inline void dit4_block(std::vector<ull> &vec, int j, int m, const ull *itwp, const ull *itwc, const ull *itwic) {
	for(int k = 0; k < m; k++)
		dit4_butterfly<M>(vec, j + k, m, itwp[k], itwc[k], itwic[k]);
}
template<class M>
static inline void dit4_block_run(std::vector<ull> &vec, int j, int m, ull iW, ull iW2, ull IMI) {
	ull p = M::R, q = M::R, r = IMI;
	for(int k = 0; k < m; k++) {
		dit4_butterfly<M>(vec, j + k, m, p, q, r);
		M::mulmod(p, iW2), M::mulmod(q, iW), M::mulmod(r, iW);
	}
}

static inline void dif2_multi(BufRef *buf, int n, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < LOG2_OMP_PIVOT) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif2_block<M2>(*buf[a].v, j, i, M2::tw(cn));
				else dif2_block<M1>(*buf[a].v, j, i, M1::tw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif2_block_run<M2>(*buf[a].v, j, i, M2::wn(cn));
				else dif2_block_run<M1>(*buf[a].v, j, i, M1::wn(cn));
			}
	}
}
static inline void dif4_multi(BufRef *buf, int n, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < LOG2_OMP_PIVOT) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif4_block<M2>(*buf[a].v, j, m, M2::tw(cn), M2::tw(cn - 1), M2::twim(cn));
				else dif4_block<M1>(*buf[a].v, j, m, M1::tw(cn), M1::tw(cn - 1), M1::twim(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif4_block_run<M2>(*buf[a].v, j, m, M2::wn(cn), M2::wn(cn - 1), M2::im());
				else dif4_block_run<M1>(*buf[a].v, j, m, M1::wn(cn), M1::wn(cn - 1), M1::im());
			}
	}
}
static inline void dit2_multi(BufRef *buf, int n, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < LOG2_OMP_PIVOT) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit2_block<M2>(*buf[a].v, j, i, M2::itw(cn));
				else dit2_block<M1>(*buf[a].v, j, i, M1::itw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit2_block_run<M2>(*buf[a].v, j, i, M2::iwn(cn));
				else dit2_block_run<M1>(*buf[a].v, j, i, M1::iwn(cn));
			}
	}
}
static inline void dit4_multi(BufRef *buf, int n, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < LOG2_OMP_PIVOT) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit4_block<M2>(*buf[a].v, j, m, M2::itw(cn - 1), M2::itw(cn), M2::itwim(cn));
				else dit4_block<M1>(*buf[a].v, j, m, M1::itw(cn - 1), M1::itw(cn), M1::itwim(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit4_block_run<M2>(*buf[a].v, j, m, M2::iwn(cn), M2::iwn(cn - 1), M2::iminv());
				else dit4_block_run<M1>(*buf[a].v, j, m, M1::iwn(cn), M1::iwn(cn - 1), M1::iminv());
			}
	}
}
static inline void dif_stages_multi(BufRef *buf, int n, int lim, int width) {
	int cn = width;
	if(width & 1) { dif2_multi(buf, n, lim, cn); cn--; }
	for(; cn >= 2; cn -= 2) dif4_multi(buf, n, lim, cn);
}
static inline void dit_stages_multi(BufRef *buf, int n, int lim, int width) {
	int cmax = (width & 1) ? width - 1 : width;
	for(int cn = 2; cn <= cmax; cn += 2) dit4_multi(buf, n, lim, cn);
	if(width & 1) dit2_multi(buf, n, lim, width);
	#pragma omp for
	for(int i = 0; i < lim; i++)
		for(int a = 0; a < n; a++) {
			std::vector<ull> &vec = *buf[a].v;
			if(buf[a].m2) M2::mulmod(vec[i], M2::invlimf(width));
			else M1::mulmod(vec[i], M1::invlimf(width));
		}
}

// ---- multiplication -------------------------------------------------------
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
		x.w.resize(lim), y.w.resize(lim), tmp1 = x.w, tmp2 = y.w;
		auto &v1 = x.w, &v2 = y.w, &v3 = tmp1, &v4 = tmp2;
#ifdef SPEED
		SpeedRec &sr = speed_rec[lim];
		++sr.calls;
#endif
		bool omp = lim >= OMP_PIVOT;
#ifdef SPEED
		std::chrono::steady_clock::time_point tprev = std::chrono::steady_clock::now();
#endif
		if(omp) {
			// One parallel region for the whole multiplication.  The M1 and M2
			// pipelines -- and the two arrays of each -- are merged into shared
			// worksharing loops so that the narrow early-DIF / late-DIT stages
			// keep all threads busy.
			BufRef dbuf[4] = {{&v3, false}, {&v4, false}, {&v1, true}, {&v2, true}};
			BufRef ditbuf[2] = {{&v3, false}, {&v1, true}};
			#pragma omp parallel
			{
				#pragma omp for
				for(int i = 0; i < lim; i++) {
					to_mont1(v3[i]), to_mont1(v4[i]);
					to_mont2(v1[i]), to_mont2(v2[i]);
				}
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(to_mont); }
#endif
				dif_stages_multi(dbuf, 4, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dif); }
#endif
				#pragma omp for
				for(int i = 0; i < lim; i++) {
					mulmod1(v3[i], v4[i]);
					mulmod2(v1[i], v2[i]);
				}
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(pointwise); }
#endif
				dit_stages_multi(ditbuf, 2, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dit); }
#endif
			}
		} else {
			SPEED_TICK(a0);
			for(int i = 0; i < lim; i++) to_mont1(v3[i]), to_mont1(v4[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			dif_ser<M1>(v3, lim, width);
			dif_ser<M1>(v4, lim, width);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) mulmod1(v3[i], v4[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			SPEED_TICK(a5);
			dit_ser<M1>(v3, lim, width);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);

			SPEED_TICK(b0);
			for(int i = 0; i < lim; i++) to_mont2(v1[i]), to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			dif_ser<M2>(v1, lim, width);
			dif_ser<M2>(v2, lim, width);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) mulmod2(v1[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			SPEED_TICK(b5);
			dit_ser<M2>(v1, lim, width);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}

		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		// (the de-Montgomery conversion is already folded into the DIT scaling)
		u128 last = 0;
		for(int i = 0; i < lim; i++) {
			ull t = (v1[i] >= v3[i] ? v1[i] - v3[i] : v1[i] + MOD2 - v3[i]);
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
		bool omp = lim >= OMP_PIVOT;
#ifdef SPEED
		std::chrono::steady_clock::time_point tprev = std::chrono::steady_clock::now();
#endif
		if(omp) {
			BufRef dbuf[2] = {{&v1, false}, {&v2, true}};
			#pragma omp parallel
			{
				#pragma omp for
				for(int i = 0; i < lim; i++) {
					to_mont1(v1[i]);
					to_mont2(v2[i]);
				}
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(to_mont); }
#endif
				dif_stages_multi(dbuf, 2, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dif); }
#endif
				#pragma omp for
				for(int i = 0; i < lim; i++) {
					mulmod1(v1[i], v1[i]);
					mulmod2(v2[i], v2[i]);
				}
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(pointwise); }
#endif
				dit_stages_multi(dbuf, 2, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dit); }
#endif
			}
		} else {
			SPEED_TICK(a0);
			for(int i = 0; i < lim; i++) to_mont1(v1[i]);
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			dif_ser<M1>(v1, lim, width);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) mulmod1(v1[i], v1[i]);
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			SPEED_TICK(a5);
			dit_ser<M1>(v1, lim, width);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);

			SPEED_TICK(b0);
			for(int i = 0; i < lim; i++) to_mont2(v2[i]);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			dif_ser<M2>(v2, lim, width);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) mulmod2(v2[i], v2[i]);
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			SPEED_TICK(b5);
			dit_ser<M2>(v2, lim, width);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}

		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		// (the de-Montgomery conversion is already folded into the DIT scaling)
		u128 last = 0;
		for(int i = 0; i < lim; i++) {
			ull t = (v2[i] >= v1[i] ? v2[i] - v1[i] : v2[i] + MOD2 - v1[i]);
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

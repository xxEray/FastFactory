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
// Twiddle rows: row w holds W_w^i for i < 2^(w-1).  Three quantities needed
// by a radix-4 butterfly at stage cn all read this single row:
//   p = W^k        = tw[k]
//   q = W^{2k}     = tw[2k]           (W_{cn-1} = W_cn^2)
//   r = IM·W^k     = tw[k + 2^(cn-2)] (IM = W^{2^(cn-2)})
// so one row replaces the tw/twp/twim trio per stage and the butterflies
// stream one table instead of three.  Stages cn < TW_MAX use these rows; the
// handful of larger stages advance their twiddles by multiplication.
// Total: 4 rows-sets x 2^21 entries x 8 B = 64 MB.
constexpr int TW_MAX = 22;
std::vector<ull> tw1[TW_MAX], tw2[TW_MAX];
std::vector<ull> itw1[TW_MAX], itw2[TW_MAX];
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
	for(int w = 0; w < TW_MAX; w++) {
		int cnt = w == 0 ? 1 : (1 << (w - 1)); // rows cover i < 2^(w-1)
		tw1[w].assign(cnt, 0), itw1[w].assign(cnt, 0);
		tw2[w].assign(cnt, 0), itw2[w].assign(cnt, 0);
		tw1[w][0] = R_MOD1, itw1[w][0] = R_MOD1;
		tw2[w][0] = R_MOD2, itw2[w][0] = R_MOD2;
		for(int i = 1; i < cnt; i++) {
			tw1[w][i] = tw1[w][i - 1], mulmod1(tw1[w][i], wn1[w]);
			itw1[w][i] = itw1[w][i - 1], mulmod1(itw1[w][i], iwn1[w]);
			tw2[w][i] = tw2[w][i - 1], mulmod2(tw2[w][i], wn2[w]);
			itw2[w][i] = itw2[w][i - 1], mulmod2(itw2[w][i], iwn2[w]);
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
	static inline ull *tw(int c) { return tw1[c].data(); }
	static inline ull *itw(int c) { return itw1[c].data(); }
	static inline ull invlim(int w) { return invlim1[w]; }
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
	static inline ull *tw(int c) { return tw2[c].data(); }
	static inline ull *itw(int c) { return itw2[c].data(); }
	static inline ull invlim(int w) { return invlim2[w]; }
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
static inline void dif4_butterfly(ull * __restrict vec, int o, int m, ull p, ull q, ull r) {
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
static inline void dit4_butterfly(ull * __restrict vec, int o, int m, ull p, ull q, ull r) {
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
void dif2_ser(ull * __restrict vec, int lim, int cn) {
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
void dif4_ser(ull * __restrict vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	const ull *tw = M::tw(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < m; k++)
			dif4_butterfly<M>(vec, j + k, m, tw[k], tw[k << 1], tw[k + m]);
}

template<class M>
void dit2_ser(ull * __restrict vec, int lim, int cn) {
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
void dit4_ser(ull * __restrict vec, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	const ull *itw = M::itw(cn);
	for(int j = 0; j < lim; j += step)
		for(int k = 0; k < m; k++)
			dit4_butterfly<M>(vec, j + k, m, itw[k << 1], itw[k], itw[k + m]);
}

template<class M>
void dif_ser(ull * __restrict vec, int lim, int width) {
	int cn = width;
	if(width & 1) { dif2_ser<M>(vec, lim, cn); cn--; }
	for(; cn >= 2; cn -= 2) dif4_ser<M>(vec, lim, cn);
}
template<class M>
void dit_ser(ull * __restrict vec, int lim, int width) {
	int cmax = (width & 1) ? width - 1 : width;
	for(int cn = 2; cn <= cmax; cn += 2) dit4_ser<M>(vec, lim, cn);
	if(width & 1) dit2_ser<M>(vec, lim, width);
}

// ---- multi-array transforms ------------------------------------------------
// Same stages as above, but a single worksharing loop covers several arrays at
// once (with `collapse(2)` over {array, block}).  This widens the parallelism
// of the early DIF / late DIT stages, which have only a handful of blocks each
// and would otherwise leave most threads idle.
// m2 selects the modulus; src, when non-null, is the second operand of the
// pointwise product that is fused into the DIT's tiled pass (src == v gives
// the in-place square).
struct BufRef { std::vector<ull> *v; bool m2; std::vector<ull> *src = nullptr; };

template<class M>
static inline void dif2_block(ull * __restrict vec, int j, int i, const ull * __restrict twc) {
	for(int k = 0; k < i; k++) {
		ull x = vec[j + k], y = vec[j + i + k];
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
		M::mulmod(vec[j + i + k], twc[k]);
	}
}
template<class M>
static inline void dif2_block_run(ull * __restrict vec, int j, int i, ull wn) {
	ull w = M::R;
	for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
		ull x = vec[j + k], y = vec[j + i + k];
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
		M::mulmod(vec[j + i + k], w);
	}
}
template<class M>
static inline void dif4_block(ull * __restrict vec, int j, int m, const ull * __restrict tw) {
	// One row serves all three twiddles: W^k = tw[k], W^{2k} = tw[2k]
	// (W_{cn-1} = W_cn^2) and IM·W^k = W^{k+2^(cn-2)} = tw[k+m].
	for(int k = 0; k < m; k++)
		dif4_butterfly<M>(vec, j + k, m, tw[k], tw[k << 1], tw[k + m]);
}
template<class M>
static inline void dif4_block_run(ull * __restrict vec, int j, int m, ull W, ull W2, ull IM) {
	ull p = M::R, q = M::R, r = IM;
	for(int k = 0; k < m; k++) {
		dif4_butterfly<M>(vec, j + k, m, p, q, r);
		M::mulmod(p, W), M::mulmod(q, W2), M::mulmod(r, W);
	}
}
template<class M>
static inline void dit2_block(ull * __restrict vec, int j, int i, const ull * __restrict itwc) {
	for(int k = 0; k < i; k++) {
		ull x = vec[j + k], y = vec[j + i + k];
		M::mulmod(y, itwc[k]);
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
	}
}
template<class M>
static inline void dit2_block_run(ull * __restrict vec, int j, int i, ull wn) {
	ull w = M::R;
	for(int k = 0; k < i; k++, M::mulmod(w, wn)) {
		ull x = vec[j + k], y = vec[j + i + k];
		M::mulmod(y, w);
		vec[j + k] = M::trim(x + y);
		vec[j + i + k] = (x >= y ? x - y : x + M::MOD - y);
	}
}
template<class M>
static inline void dit4_block(ull * __restrict vec, int j, int m, const ull * __restrict itw) {
	// Mirror of dif4_block: W_{cn-1}^{-k} = itw[2k], W^{-k} = itw[k],
	// IM^{-1}·W^{-k} = W^{-(k+m)} = itw[k+m].
	for(int k = 0; k < m; k++)
		dit4_butterfly<M>(vec, j + k, m, itw[k << 1], itw[k], itw[k + m]);
}
template<class M>
static inline void dit4_block_run(ull * __restrict vec, int j, int m, ull iW, ull iW2, ull IMI) {
	ull p = M::R, q = M::R, r = IMI;
	for(int k = 0; k < m; k++) {
		dit4_butterfly<M>(vec, j + k, m, p, q, r);
		M::mulmod(p, iW2), M::mulmod(q, iW), M::mulmod(r, iW);
	}
}

static inline void dif2_multi(BufRef *buf, int n, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < TW_MAX) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif2_block<M2>(buf[a].v->data(), j, i, M2::tw(cn));
				else dif2_block<M1>(buf[a].v->data(), j, i, M1::tw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif2_block_run<M2>(buf[a].v->data(), j, i, M2::wn(cn));
				else dif2_block_run<M1>(buf[a].v->data(), j, i, M1::wn(cn));
			}
	}
}
static inline void dif4_multi(BufRef *buf, int n, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < TW_MAX) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif4_block<M2>(buf[a].v->data(), j, m, M2::tw(cn));
				else dif4_block<M1>(buf[a].v->data(), j, m, M1::tw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dif4_block_run<M2>(buf[a].v->data(), j, m, M2::wn(cn), M2::wn(cn - 1), M2::im());
				else dif4_block_run<M1>(buf[a].v->data(), j, m, M1::wn(cn), M1::wn(cn - 1), M1::im());
			}
	}
}
static inline void dit2_multi(BufRef *buf, int n, int lim, int cn) {
	int i = 1 << (cn - 1), step = i << 1;
	if(cn < TW_MAX) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit2_block<M2>(buf[a].v->data(), j, i, M2::itw(cn));
				else dit2_block<M1>(buf[a].v->data(), j, i, M1::itw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit2_block_run<M2>(buf[a].v->data(), j, i, M2::iwn(cn));
				else dit2_block_run<M1>(buf[a].v->data(), j, i, M1::iwn(cn));
			}
	}
}
static inline void dit4_multi(BufRef *buf, int n, int lim, int cn) {
	int m = 1 << (cn - 2), step = m << 2;
	if(cn < TW_MAX) {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit4_block<M2>(buf[a].v->data(), j, m, M2::itw(cn));
				else dit4_block<M1>(buf[a].v->data(), j, m, M1::itw(cn));
			}
	} else {
		#pragma omp for collapse(2)
		for(int a = 0; a < n; a++)
			for(int j = 0; j < lim; j += step) {
				if(buf[a].m2) dit4_block_run<M2>(buf[a].v->data(), j, m, M2::iwn(cn), M2::iwn(cn - 1), M2::iminv());
				else dit4_block_run<M1>(buf[a].v->data(), j, m, M1::iwn(cn), M1::iwn(cn - 1), M1::iminv());
			}
	}
}
// ---- cache-tiled small stages --------------------------------------------
// Stages cn <= TILE_CN only touch indices inside each 2^cn block, and blocks
// of that size nest inside a TILE-sized chunk.  Running all of those stages
// per chunk (in their usual relative order) means the chunk is fetched from
// memory once instead of once per stage: for a 2^26 transform this turns 13
// full-array passes into 6 (5 big-stride stages + 1 tiled pass).
constexpr int TILE_CN = 14;
constexpr int TILE = 1 << TILE_CN;

static inline void dif4_tiled_multi(BufRef *buf, int n, int lim, int cn) {
	#pragma omp for collapse(2)
	for(int a = 0; a < n; a++)
		for(int t = 0; t < lim; t += TILE) {
			std::vector<ull> &vec = *buf[a].v;
			int tend = t + TILE > lim ? lim : t + TILE;
			for(int c = cn; c >= 2; c -= 2) {
				int m = 1 << (c - 2), step = m << 2;
				for(int j = t; j < tend; j += step) {
					if(buf[a].m2) dif4_block<M2>(vec.data(), j, m, M2::tw(c));
					else dif4_block<M1>(vec.data(), j, m, M1::tw(c));
				}
			}
		}
}
static inline void dit4_tiled_multi(BufRef *buf, int n, int lim, int cmax) {
	#pragma omp for collapse(2)
	for(int a = 0; a < n; a++)
		for(int t = 0; t < lim; t += TILE) {
			std::vector<ull> &vec = *buf[a].v;
			int tend = t + TILE > lim ? lim : t + TILE;
			if(buf[a].src) {
				// Pointwise product of the two transformed operands, folded
				// into this pass so the arrays are streamed once less.  The
				// second operand is converted to Montgomery form on the fly:
				// mulmod(x, mont(y)) = x*y stays plain.
				std::vector<ull> &sv = *buf[a].src;
				if(buf[a].m2) {
					for(int i = t; i < tend; i++) {
						ull b = sv[i]; to_mont2(b); mulmod2(vec[i], b);
					}
				} else {
					for(int i = t; i < tend; i++) {
						ull b = sv[i]; to_mont1(b); mulmod1(vec[i], b);
					}
				}
			}
			for(int c = 2; c <= cmax; c += 2) {
				int m = 1 << (c - 2), step = m << 2;
				for(int j = t; j < tend; j += step) {
					if(buf[a].m2) dit4_block<M2>(vec.data(), j, m, M2::itw(c));
					else dit4_block<M1>(vec.data(), j, m, M1::itw(c));
				}
			}
		}
}
static inline void dif_stages_multi(BufRef *buf, int n, int lim, int width) {
	int cn = width;
	if(width & 1) { dif2_multi(buf, n, lim, cn); cn--; }
	for(; cn > TILE_CN; cn -= 2) dif4_multi(buf, n, lim, cn);
	if(cn >= 2) dif4_tiled_multi(buf, n, lim, cn);
}
static inline void dit_stages_multi(BufRef *buf, int n, int lim, int width) {
	int cmax = (width & 1) ? width - 1 : width;
	if(cmax >= 2) {
		int tmax = cmax < TILE_CN ? cmax : TILE_CN;
		if(tmax >= 2) dit4_tiled_multi(buf, n, lim, tmax);
		for(int cn = tmax + 2; cn <= cmax; cn += 2) dit4_multi(buf, n, lim, cn);
	}
	if(width & 1) dit2_multi(buf, n, lim, width);
}

// ---- multiplication -------------------------------------------------------
// Reconstruct the base-B digits of an exact product from its residues modulo
// MOD1 (r1) and MOD2 (r2), folding in the final division by 2^width.  The
// exact convolution coefficient at position i is < B^(lim-i), so splitting it
// as r + q1*B + q2*B^2 (q2 < 2^31) keeps the normalization carry <= 2: the
// heavy u128 division work runs in parallel and the serial part is a cheap
// add/compare chain.  `out` may alias r1 or r2; q1s/q2s may alias either
// residue array (each element is read before it is overwritten).
static void crt_finalize(ull *out, ull *r1, ull *r2, ull *q1s, ull *q2s,
	int lim, ull inv1, ull inv2, bool par) {
	#pragma omp parallel for schedule(static) if(par)
	for(int i = 0; i < lim; i++) {
		ull a = r1[i], b = r2[i];
		mulmod1(a, inv1), mulmod2(b, inv2);
		ull t = (b >= a ? b - a : b + MOD2 - a);
		mulmod2(t, INV_MOD1_MOD2_R);
		u128 val = static_cast<u128>(a) + static_cast<u128>(MOD1) * t;
		u128 q = val / B;
		q1s[i] = static_cast<ull>(q % B);
		q2s[i] = static_cast<ull>(q / B);
		out[i] = static_cast<ull>(val - q * B);
	}
	ull carry = 0, q1m = 0, q2m1 = 0, q2m2 = 0;
	for(int i = 0; i < lim; i++) {
		ull s = out[i] + q1m + q2m2 + carry;
		carry = 0;
		if(s >= B) { s -= B, carry = 1; }
		if(s >= B) { s -= B, carry = 2; }
		out[i] = s;
		q2m2 = q2m1, q2m1 = q2s[i], q1m = q1s[i];
	}
}

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
			//
			// Values stay in plain (non-Montgomery) form throughout: mulmod of
			// a plain value by a Montgomery-form twiddle yields W*x, i.e. the
			// Montgomery factors cancel.  The conversion of the second operand
			// is folded into the pointwise product, which itself is folded
			// into the DIT's tiled pass, and the de-Montgomery into the CRT
			// pass: no separate array sweeps remain between DIF and DIT.
			BufRef dbuf[4] = {{&v3, false}, {&v4, false}, {&v1, true}, {&v2, true}};
			BufRef ditbuf[2] = {{&v3, false, &v4}, {&v1, true, &v2}};
			#pragma omp parallel
			{
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(to_mont); }
#endif
				dif_stages_multi(dbuf, 4, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dif); }
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
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			dif_ser<M1>(v3.data(), lim, width);
			dif_ser<M1>(v4.data(), lim, width);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) { ull a = v4[i]; to_mont1(a); mulmod1(v3[i], a); }
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			SPEED_TICK(a5);
			dit_ser<M1>(v3.data(), lim, width);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);

			SPEED_TICK(b0);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			dif_ser<M2>(v1.data(), lim, width);
			dif_ser<M2>(v2.data(), lim, width);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) { ull b = v2[i]; to_mont2(b); mulmod2(v1[i], b); }
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			SPEED_TICK(b5);
			dit_ser<M2>(v1.data(), lim, width);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}

		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).  The u128 work runs
		// in parallel and the normalization carries are bounded by 2, so the
		// serial part is a cheap add/compare chain.  v2/v4 are dead after the
		// pointwise pass and carry the digit streams; the final division by
		// 2^width (formerly a separate pass after DIT) is folded in here.
		crt_finalize(v1.data(), v3.data(), v1.data(), v2.data(), v4.data(),
			lim, M1::invlim(width), M2::invlim(width), omp);
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
			// src == v: the pointwise square is fused into the DIT's tiled
			// pass (a Montgomery copy of each element supplies the R factor).
			BufRef dbuf[2] = {{&v1, false, &v1}, {&v2, true, &v2}};
			#pragma omp parallel
			{
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(to_mont); }
#endif
				dif_stages_multi(dbuf, 2, lim, width);
#ifdef SPEED
				#pragma omp master
				{ SPEED_MARK(dif); }
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
			SPEED_TICK(a1);
			SPEED_ADD(to_mont, a0, a1);
			dif_ser<M1>(v1.data(), lim, width);
			SPEED_TICK(a2);
			SPEED_ADD(dif, a1, a2);
			SPEED_TICK(a3);
			for(int i = 0; i < lim; i++) { ull a = v1[i]; to_mont1(a); mulmod1(v1[i], a); }
			SPEED_TICK(a4);
			SPEED_ADD(pointwise, a3, a4);
			SPEED_TICK(a5);
			dit_ser<M1>(v1.data(), lim, width);
			SPEED_TICK(a6);
			SPEED_ADD(dit, a5, a6);

			SPEED_TICK(b0);
			SPEED_TICK(b1);
			SPEED_ADD(to_mont, b0, b1);
			dif_ser<M2>(v2.data(), lim, width);
			SPEED_TICK(b2);
			SPEED_ADD(dif, b1, b2);
			SPEED_TICK(b3);
			for(int i = 0; i < lim; i++) { ull b = v2[i]; to_mont2(b); mulmod2(v2[i], b); }
			SPEED_TICK(b4);
			SPEED_ADD(pointwise, b3, b4);
			SPEED_TICK(b5);
			dit_ser<M2>(v2.data(), lim, width);
			SPEED_TICK(b6);
			SPEED_ADD(dit, b5, b6);
		}

		SPEED_TICK(t8);
		// CRT: reconstruct c in [0, MOD1*MOD2) from c mod MOD1 and c mod MOD2.
		// x = r1 + MOD1 * (((r2-r1) * inv(MOD1) mod MOD2)).
		// v2 is consumed here and carries the q1 digit stream; tmp2 (unused by
		// this routine until now) carries q2.
		tmp2.resize(lim);
		crt_finalize(v1.data(), v1.data(), v2.data(), v2.data(), tmp2.data(),
			lim, M1::invlim(width), M2::invlim(width), omp);
		SPEED_TICK(t9);
		SPEED_ADD(crt, t8, t9);
		x.pop_zero();
	}
}

#endif

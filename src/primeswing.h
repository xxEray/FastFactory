#ifndef _PRIMESWING_
#define _PRIMESWING_

#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include "common.h"
#include "bigint.h"

// Odd-only bit-packed sieve.  Bit (p >> 1) stores the primality of the odd
// number p; 2 is handled specially.  Memory drops from n bytes (~100 MB at
// n = 1e8) to n/16 bytes (~6 MB), which keeps the working set cache-friendly
// and leaves headroom under the 3 GB limit.
std::vector<ull> prm_bits;

static inline bool isprm(int p) {
	if(p < 2) return false;
	if(p == 2) return true;
	if(!(p & 1)) return false;
	return (prm_bits[p >> 7] >> ((p >> 1) & 63)) & 1;
}

void sieve(int n) {
	SPEED_TICK(s0);
	const size_t half = (size_t)n >> 1; // highest index = n >> 1
	prm_bits.assign(half / 64 + 2, ~0ULL);
	prm_bits[0] &= ~1ULL; // 1 is not prime
	for(int q = 3; (ull)q * q <= (ull)n; q += 2) {
		if(!isprm(q)) continue;
		for(size_t j = ((size_t)q * q) >> 1; j <= half; j += q)
			prm_bits[j >> 6] &= ~(1ULL << (j & 63));
	}
	SPEED_TICK(s1);
	SPEED_ADD_TOTAL(sieve, s0, s1);
}

// Levelized product tree.
//
// The original implementation was a plain recursion:
//
//     val = product_tree(l, mid); mul_eq(val, product_tree(mid + 1, r));
//
// i.e. the two child products were built strictly one after the other, so the
// whole lower part of the tree ran single-threaded (lim < OMP_PIVOT => the NTT
// runs serially).
//
// Here we build LEAF_PIVOT-sized leaf products in parallel, then reduce them
// pairwise.  A round is parallelized across pairs while it still has at least
// TREE_PAR_MIN pairs: each pair then runs its NTT serially (no nested region),
// which keeps all threads busy even though a single transform only scales
// ~1.7x internally.  Once the pair count drops below TREE_PAR_MIN the round is
// run serially so the last few big multiplications can use the whole team.
//
// Memory: at a level with `nxt` pairs the operands total ~|result| limbs, so
// the concurrent scratch is <= 32 * |result| bytes -- i.e. no more than what a
// single top-level multiplication already needs.  The product itself is
// unchanged.
constexpr int TREE_PAR_MIN = 2;

BigInt product_tree(const std::vector<unsigned> &vec, int l, int r) {
	if(l > r) return 1u;
	int n = r - l + 1;
	if(n <= LEAF_PIVOT) {
		BigInt ret(vec[l]);
		for(int i = l + 1; i <= r; i++) ret *= vec[i];
		return ret;
	}
	int nblk = (n + LEAF_PIVOT - 1) / LEAF_PIVOT;
	std::vector<BigInt> cur(nblk);
	#pragma omp parallel for schedule(dynamic)
	for(int b = 0; b < nblk; b++) {
		int lo = l + b * LEAF_PIVOT;
		int hi = std::min(r, lo + LEAF_PIVOT - 1);
		BigInt t(vec[lo]);
		for(int i = lo + 1; i <= hi; i++) t *= vec[i];
		cur[b] = std::move(t);
	}
	while(cur.size() > 1) {
		int m = (int)cur.size(), nxt = (m + 1) / 2;
		std::vector<BigInt> nxtv(nxt);
		if(nxt >= TREE_PAR_MIN) {
			#pragma omp parallel for schedule(dynamic)
			for(int i = 0; i < nxt; i++) {
				if(2 * i + 1 < m) {
					BigInt t = std::move(cur[2 * i]);
					mul_eq(t, std::move(cur[2 * i + 1]));
					nxtv[i] = std::move(t);
				} else nxtv[i] = std::move(cur[2 * i]);
			}
		} else {
			for(int i = 0; i < nxt; i++) {
				if(2 * i + 1 < m) {
					BigInt t = std::move(cur[2 * i]);
					mul_eq(t, std::move(cur[2 * i + 1]));
					nxtv[i] = std::move(t);
				} else nxtv[i] = std::move(cur[2 * i]);
			}
		}
		cur.swap(nxtv);
	}
	return std::move(cur[0]);
}

BigInt swing(int n) {
	SPEED_TICK(w0);
	std::vector<unsigned> vec = {1u << __builtin_popcount(n / 2)};
	int sqn = static_cast<int>(sqrtl(n));
	while(static_cast<ull>(sqn + 1) * (sqn + 1) <= n) sqn++;
	while(static_cast<ull>(sqn) * sqn > n) sqn--;
	assert(n >= 4);
	for(int p = 3; p <= sqn; p++) if(isprm(p)) {
		int t = n;
		unsigned val = 1;
		while(t) {
			t /= p;
			if(t & 1) val *= p;
		}
		if(val > 1) vec.emplace_back(val);
	}
	for(int p = sqn + 1; p <= n / 3; p++) if(isprm(p))
		if(n / p & 1) vec.emplace_back(p);
	for(int p = n / 2 + 1; p <= n; p++) if(isprm(p)) vec.emplace_back(p);
	SPEED_TICK(w1);
	SPEED_ADD_TOTAL(swing, w0, w1);
	BigInt ret = product_tree(vec, 0, (int)vec.size() - 1);
	return ret;
}

BigInt factorial(int n) {
	if(n <= 1) return 1u;
	if(n <= FAC_BF_PIVOT) {
		BigInt val = 1u;
		for(int i = 2; i <= n; i++) val *= static_cast<unsigned>(i);
		return val;
	}
	// n! = f(n/2)^2 * swing(n) = f(n/2) * (f(n/2) * swing(n)).  Folding the
	// swing into one factor first keeps that product short (|f|+|swing| instead
	// of 2|f|), so the first multiplication transforms at half the length; the
	// second then replaces the squaring.  In transform work this is
	// 6*N/2*log(N/2) + 6*N*log(N) instead of 4*N*log(N) + 6*N*log(N).
	BigInt val = factorial(n / 2);
	BigInt t = val;
	mul_eq(t, swing(n));
	mul_eq(val, std::move(t));
	return val;
}

#endif
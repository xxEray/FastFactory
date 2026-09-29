#ifndef _PRIMESWING_
#define _PRIMESWING_

#include <cmath>
#include <vector>
#include <functional>
#include "common.h"
#include "bigint.h"

bool *isprm;
void sieve(int n) {
	SPEED_TICK(s0);
	isprm = new bool[n + 5];
	std::fill(isprm, isprm + n + 5, true);
	isprm[0] = isprm[1] = false;
	for(int i = 2; i <= n; i++) if(isprm[i] && static_cast<ull>(i) * i <= n)
		for(int j = i * i; j <= n; j += i) isprm[j] = false;
	SPEED_TICK(s1);
	SPEED_ADD_TOTAL(sieve, s0, s1);
}

BigInt product_tree(const std::vector<unsigned> &vec, int l, int r) {
	if(l > r) return 1u;
	if(r - l + 1 <= LEAF_PIVOT) {
		BigInt ret(vec[l]);
		for(int i = l + 1; i <= r; i++) ret *= vec[i];
		return ret;
	}
	int mid = (l + r) / 2;
	BigInt val = product_tree(vec, l, mid);
	mul_eq(val, product_tree(vec, mid + 1, r));
	return val;
}

BigInt swing(int n) {
	SPEED_TICK(w0);
	std::vector<unsigned> vec = {1u << __builtin_popcount(n / 2)};
	int sqn = static_cast<int>(sqrtl(n));
	while(static_cast<ull>(sqn + 1) * (sqn + 1) <= n) sqn++;
	while(static_cast<ull>(sqn) * sqn > n) sqn--;
	assert(n >= 4);
	for(int p = 3; p <= sqn; p++) if(isprm[p]) {
		int t = n;
		unsigned val = 1;
		while(t) {
			t /= p;
			if(t & 1) val *= p;
		}
		if(val > 1) vec.emplace_back(val);
	}
	for(int p = sqn + 1; p <= n / 3; p++) if(isprm[p])
		if(n / p & 1) vec.emplace_back(p);
	for(int p = n / 2 + 1; p <= n; p++) if(isprm[p]) vec.emplace_back(p);
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
	BigInt val = factorial(n / 2);
	mul_self_eq(val);
	mul_eq(val, swing(n));
	return val;
}

#endif
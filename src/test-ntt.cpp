#include <iostream>
#include <string>

#include "bigint.h"

BigInt parse(const std::string &s) {
	BigInt r;
	r.w.clear();
	for(int i = s.size(); i > 0; i -= 14) {
		int l = std::max(0, i - 14);
		ull v = 0;
		for(int j = l; j < i; j++) v = v * 10 + (s[j] - '0');
		r.w.push_back(v);
	}
	if(r.w.empty()) r.w.push_back(0);
	r.pop_zero();
	return r;
}

std::string print(const BigInt &x) {
	std::string s;
	for(int i = x.w.size() - 1; i >= 0; i--) {
		std::string t = std::to_string(x.w[i]);
		if(i != static_cast<int>(x.w.size()) - 1)
			t = std::string(14 - t.size(), '0') + t;
		s += t;
	}
	return s;
}

int main() {
	std::string a, b;
	if(!(std::cin >> a >> b)) return 0;
	BigInt ia = parse(a), ib = parse(b);
	mul_eq(ia, std::move(ib));
	std::cout << print(ia) << "\n";
	return 0;
}

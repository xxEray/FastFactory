#ifndef _SPEED_
#define _SPEED_

#ifdef SPEED
#include <chrono>
#include <map>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#endif

struct SpeedRec {
	long long to_mont = 0, dif = 0, pointwise = 0, dit = 0, crt = 0, calls = 0;
};
std::map<int, SpeedRec> speed_rec, speed_rec_self;

struct SpeedDump {
	long long sieve = 0, swing = 0, prod = 0, output = 0;
	std::chrono::steady_clock::time_point start;
	SpeedDump() : start(std::chrono::steady_clock::now()) {}
	static void dump(FILE *fp, const char *title, const std::map<int, SpeedRec> &rec) {
		fprintf(fp, "[%s]\n", title);
		fprintf(fp, "%12s%16s%16s%16s%16s%16s%14s\n", "lim", "to_mont(ms)", "dif(ms)", "point(ms)", "dit(ms)", "crt(ms)", "calls");
		long long to_mont = 0, dif = 0, pointwise = 0, dit = 0, crt = 0, calls = 0;
		for(const auto &it : rec) {
			if(it.first == 0) fprintf(fp, "%12s", "small");
			else fprintf(fp, "%12d", it.first);
			fprintf(fp, "%16.3f%16.3f%16.3f%16.3f%16.3f%14lld\n", it.second.to_mont / 1e6,
				it.second.dif / 1e6, it.second.pointwise / 1e6, it.second.dit / 1e6,
				it.second.crt / 1e6, it.second.calls);
			to_mont += it.second.to_mont, dif += it.second.dif;
			pointwise += it.second.pointwise, dit += it.second.dit;
			crt += it.second.crt, calls += it.second.calls;
		}
		fprintf(fp, "%12s%16.3f%16.3f%16.3f%16.3f%16.3f%14lld\n", "TOTAL", to_mont / 1e6,
			dif / 1e6, pointwise / 1e6, dit / 1e6, crt / 1e6, calls);
	}
	~SpeedDump() {
		FILE *fp = fopen("speed.txt", "w");
		if(!fp) return;
		dump(fp, "mul_eq", speed_rec);
		dump(fp, "mul_self_eq", speed_rec_self);
		fprintf(fp, "sieve(ms) %.3f\n", sieve / 1e6);
		fprintf(fp, "swing(ms) %.3f\n", swing / 1e6);
		fprintf(fp, "prod(ms) %.3f\n", prod / 1e6);
		fprintf(fp, "output(ms) %.3f\n", output / 1e6);
		auto end = std::chrono::steady_clock::now();
		double total_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();
		fprintf(fp, "total(ms) %.3f\n", total_ms);
#if defined(_WIN32)
		PROCESS_MEMORY_COUNTERS pmc;
		if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			fprintf(fp, "peakmem(MB) %.3f\n", pmc.PeakWorkingSetSize / 1024.0 / 1024.0);
#endif
		fclose(fp);
	}
};
SpeedDump speed_dump;

#define SPEED_TICK(name) auto name = std::chrono::steady_clock::now()
#define SPEED_ADD(field, from, to) sr.field += \
	std::chrono::duration_cast<std::chrono::nanoseconds>((to) - (from)).count()
#define SPEED_ADD_TOTAL(field, from, to) speed_dump.field += \
	std::chrono::duration_cast<std::chrono::nanoseconds>((to) - (from)).count()
#else
#define SPEED_TICK(name) ((void)0)
#define SPEED_ADD(field, from, to) ((void)0)
#define SPEED_ADD_TOTAL(field, from, to) ((void)0)
#endif

#endif
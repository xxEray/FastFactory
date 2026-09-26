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
	long long ntt = 0, intt = 0, crt = 0, calls = 0;
};
std::map<int, SpeedRec> speed_rec;

struct SpeedDump {
	std::chrono::steady_clock::time_point start;
	SpeedDump() : start(std::chrono::steady_clock::now()) {}
	~SpeedDump() {
		FILE *fp = fopen("speed.txt", "w");
		if(!fp) return;
		fprintf(fp, "%12s%16s%16s%16s%14s\n", "lim", "ntt(ms)", "intt(ms)", "crt(ms)", "calls");
		long long ntt = 0, intt = 0, crt = 0, calls = 0;
		for(const auto &it : speed_rec) {
			if(it.first == 0) fprintf(fp, "%12s", "small");
			else fprintf(fp, "%12d", it.first);
			fprintf(fp, "%16.3f%16.3f%16.3f%14lld\n", it.second.ntt / 1e6,
				it.second.intt / 1e6, it.second.crt / 1e6, it.second.calls);
			ntt += it.second.ntt, intt += it.second.intt;
			crt += it.second.crt, calls += it.second.calls;
		}
		fprintf(fp, "%12s%16.3f%16.3f%16.3f%14lld\n", "TOTAL", ntt / 1e6, intt / 1e6, crt / 1e6, calls);
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
#else
#define SPEED_TICK(name) ((void)0)
#define SPEED_ADD(field, from, to) ((void)0)
#endif

#endif
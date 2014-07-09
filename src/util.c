#include "util.h"

/**
 * Starts the high resolution timer
 */
void start_timer(struct hrtimer_t *timer) {
#ifdef LINUX
	clock_gettime(CLOCK_REALTIME, &timer->start);
#else
	QueryPerformanceFrequency(&timer->freq);
	QueryPerformanceCounter(&timer->start);
#endif
}

/**
 * Stops the high resolution timer
 */
void stop_timer(struct hrtimer_t *timer) {
#ifdef LINUX
	clock_gettime(CLOCK_REALTIME, &timer->stop);
#else
	QueryPerformanceCounter(&timer->stop);
#endif
}

/**
 * Reads the high resolution timer in seconds
 */
double get_timer_interval(struct hrtimer_t *timer) {
#ifdef LINUX
	long dnsec = timer->stop.tv_nsec - timer->start.tv_nsec;
	int dsec = timer->stop.tv_sec - timer->start.tv_sec;

	if (dnsec < 0) {
		dnsec += 1e9;
		dsec -= 1;
	}
	
	return ((double)dsec) + dnsec / 1.e-9;
#else
	return ((double)(timer->stop.QuadPart - timer->start.QuadPart))/timer->freq.QuadPart;
#endif
}

/**
 * Finds the ceiling of the log2 of a number iteratively
 */
int cb_log2(int x) {
	int res = 0;
	int x2 = x;

	while (x2 > 1) {
		x2 >>= 1;
		res += 1;
	}

	if ((1 << res) == x)
		return res;
	return res+1;
}

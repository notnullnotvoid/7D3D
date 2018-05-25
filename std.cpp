#include <random>

float random_float() {
	static std::random_device rd;
	static std::default_random_engine re;
	static std::uniform_real_distribution<float> urd(0, 1);

	return urd(re);
}

float random_float(float max) {
	return random_float() * max;
}

float random_float(float min, float max) {
	return random_float() * (max - min) + min;
}

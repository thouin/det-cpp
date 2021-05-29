#include <cmath>

#include "det.hpp"

namespace {
	I gcd(I a, I b) {
		if (!b)
			return a;
		else {
			return gcd(b, a % b);
		}
	}
}

void Rational::simplify() {
	if (m_n==I(0)) {
		m_d = I(1);
	} else {
		bool neg = (std::signbit(m_n) != std::signbit(m_d));
		I a = gcd(std::abs(m_n), std::abs(m_d));
		m_n = std::abs(m_n) / a;
		m_d = std::abs(m_d) / a;
		if (neg) m_n = -m_n;
	}
}

Rational::Rational(const I& n, const I& d) {
	assert(d);
	I cd(gcd(n, d));
	m_n = n/cd;
	m_d = d/cd;
}

Rational::Rational(I&& n, I&& d) {
	I cd(gcd(n, d));
	n /= cd;
	d /= cd;
	m_n = std::move(n);
	m_d = std::move(d);
}

#if USE_SQRT2 or USE_CX

std::ostream& operator<<(std::ostream& out, const Sqrt2Ring val) {
	bool r = val.re();
	bool i = val.im();
	if (r && i)
		out << '(' << val.re() << ") + (" << val.im() << ")*sqrt(2)";
	else if (i)
		out << '(' << val.im() << ")*sqrt(2)";
	else
		out << val.re();

	return out;
}

#endif

#if USE_CX

std::ostream& operator<<(std::ostream& out, const CxSqrt2Ring val) {
	bool r = val.re();
	bool i = val.im();
	if (r && i)
		out << '(' << val.re() << ") + (" << val.im() << ")*I";
	else if (i)
		out << '(' << val.im() << ")*I";
	else
		out << val.re();

	return out;
}

#endif

std::ostream& operator<<(std::ostream& out, const Rational& r) {
	if (r.denominator() == R(1))
		out << r.numerator();
	else
		out << r.numerator() << '/' << r.denominator();
	return out;
}

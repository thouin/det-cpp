#ifndef DET_HPP
#define DET_HPP

#define USE_INT 1
#define USE_SQRT2 0
#define USE_CX 0

#define MAX_SIZE 6
#define PARALLEL_THRESHOLD 4
#define PRINT_SKIP 511


#if USE_INT

# if USE_SQRT2 or USE_CX
#  error "Must only choose one ring"
#endif

#elif USE_SQRT2

# if USE_CX or USE_INT
#  error "Must only choose one ring"
#endif

#elif USE_CX

# if USE_SQRT2 or USE_INT
#  error "Must only choose one ring"
#endif

#else

# error "Must choose a ring"

#endif



#define MAX_MATRIX_SIZE (MAX_SIZE*(MAX_SIZE+1)/2)

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <initializer_list>
#include <array>
#include <list>
#include <iostream>
#include <cassert>

using I = int64_t;

class Rational {
	I m_n;
	I m_d;

	void simplify();

public:
	Rational() : m_n(I(0)), m_d(I(1)) {}
	Rational(const I& n) : m_n(n), m_d(I(1)) {}
	Rational(I&& n) : m_n(std::move(n)), m_d(I(1)) {}
	Rational(const I& n, const I& d);
	Rational(I&& n, I&& d);
	Rational operator-() const {return Rational(-m_n, m_d);}
	Rational operator+(const Rational& other) const {
		if (other && *this) {
			Rational ret(m_n*other.m_d + m_d*other.m_n, m_d*other.m_d);
			ret.simplify();
			return ret;
		} else if (!other) {
			return *this;
		} else {
			return other;
		}
	}
	Rational operator-(const Rational& other) const {
		if (other && *this) {
			Rational ret(m_n*other.m_d - m_d*other.m_n, m_d*other.m_d);
			ret.simplify();
			return ret;
		} else if (!other) {
			return *this;
		} else {
			return -other;
		}
	}
	Rational operator*(const Rational& other) const {
		Rational ret(m_n*other.m_n, m_d*other.m_d);
		ret.simplify();
		return ret;
	}
	Rational operator/(const Rational& other) const {
		assert(other.m_n != 0);
		Rational ret(m_n*other.m_d, m_d*other.m_n);
		ret.simplify();
		return ret;
	}
	Rational& operator+=(const Rational& other) {Rational r=(*this)+other; *this = std::move(r); return *this;}
	Rational& operator-=(const Rational& other) {Rational r=(*this)-other; *this = std::move(r); return *this;}
	Rational& operator*=(const Rational& other) {m_n *= other.m_n; m_d *= other.m_d; simplify(); return *this;}
	Rational& operator/=(const Rational& other) {assert(other.m_n!=I(0)); m_n *= other.m_d; m_d *= other.m_n; simplify(); return *this;}
	operator bool() const {return static_cast<bool>(m_n);}
	const I& numerator() const {return m_n;}
	const I& denominator() const {return m_d;}
};

#if USE_SQRT2 or USE_CX

using J = Rational;

class Sqrt2Ring {
	J m_r;
	J m_i;

public:
	Sqrt2Ring() =default;
	Sqrt2Ring(const J& r, const J& i) : m_r(r), m_i(i) {}
	explicit Sqrt2Ring(const J& r) : m_r(r), m_i(J()) {}
	J& re() {return m_r;}
	J& im() {return m_i;}
	const J& re() const {return m_r;}
	const J& im() const {return m_i;}
	void conj() {m_i = -m_i;}
	Sqrt2Ring getConj() const {return Sqrt2Ring(m_r, -m_i);}
	Sqrt2Ring operator+(const Sqrt2Ring& other) const {return Sqrt2Ring(m_r+other.m_r, m_i+other.m_i);}
	Sqrt2Ring operator-(const Sqrt2Ring& other) const {return Sqrt2Ring(m_r-other.m_r, m_i-other.m_i);}
	Sqrt2Ring operator-() const {return Sqrt2Ring(-m_r, -m_i);}
	Sqrt2Ring& operator+=(const Sqrt2Ring& other) {m_r += other.m_r; m_i += other.m_i; return *this;}
	Sqrt2Ring& operator-=(const Sqrt2Ring& other) {m_r -= other.m_r; m_i -= other.m_i; return *this;}
	Sqrt2Ring operator*(const Sqrt2Ring& other) const {
		return Sqrt2Ring(m_r*other.m_r + m_i*other.m_i*J(2), m_i*other.m_r + m_r*other.m_i);
	}
	Sqrt2Ring& operator*=(const Sqrt2Ring& other) {
		Sqrt2Ring s = (*this)*other;
		*this = std::move(s);
		return *this;
	}
	Sqrt2Ring& operator/=(const Sqrt2Ring& other) {
		Sqrt2Ring s = (*this)/other;
		*this = std::move(s);
		return *this;
	}
	Sqrt2Ring operator/(const Sqrt2Ring& other) const {
		return (operator*(other.getConj())).operator/(J(other.re()*other.re() - J(2)*other.im()*other.im()));
	}
	Sqrt2Ring& operator*=(const J& i) {m_r *= i; m_i *= i; return *this;}
	Sqrt2Ring& operator/=(const J& i) {m_r /= i; m_i /= i; return *this;}
	Sqrt2Ring operator*(const J& i) const {Sqrt2Ring ret(*this); return ret *= i;}
	Sqrt2Ring operator/(const J& i) const {Sqrt2Ring ret(*this); return ret /= i;}
	operator bool() const {return static_cast<bool>(m_r) || static_cast<bool>(m_i);}
};

#endif

#if USE_CX

class CxSqrt2Ring {
	Sqrt2Ring m_r;
	Sqrt2Ring m_i;

public:
	CxSqrt2Ring() =default;
	CxSqrt2Ring(const Sqrt2Ring& r, const Sqrt2Ring& i) : m_r(r), m_i(i) {}
	explicit CxSqrt2Ring(const Sqrt2Ring& r) : m_r(r), m_i(Sqrt2Ring()) {}
	explicit CxSqrt2Ring(const Rational& r) : m_r(r), m_i(0) {}
	CxSqrt2Ring(const Rational& r, const Sqrt2Ring& i) : m_r(r), m_i(i) {}
	CxSqrt2Ring(const Sqrt2Ring& r, const Rational& i) : m_r(r), m_i(i) {}
	CxSqrt2Ring(const Rational& r, const Rational& i) : m_r(r), m_i(i) {}
	Sqrt2Ring& re() {return m_r;}
	Sqrt2Ring& im() {return m_i;}
	const Sqrt2Ring& re() const {return m_r;}
	const Sqrt2Ring& im() const {return m_i;}
	void conj() {m_i = -m_i;}
	CxSqrt2Ring getConj() const {return CxSqrt2Ring(m_r, -m_i);}
	CxSqrt2Ring operator+(const CxSqrt2Ring& other) const {return CxSqrt2Ring(m_r+other.m_r, m_i+other.m_i);}
	CxSqrt2Ring operator-(const CxSqrt2Ring& other) const {return CxSqrt2Ring(m_r-other.m_r, m_i-other.m_i);}
	CxSqrt2Ring operator-() const {return CxSqrt2Ring(-m_r, -m_i);}
	CxSqrt2Ring& operator+=(const CxSqrt2Ring& other) {m_r += other.m_r; m_i += other.m_i; return *this;}
	CxSqrt2Ring& operator-=(const CxSqrt2Ring& other) {m_r -= other.m_r; m_i -= other.m_i; return *this;}
	CxSqrt2Ring operator*(const CxSqrt2Ring& other) const {
		return CxSqrt2Ring(m_r*other.m_r - m_i*other.m_i, m_i*other.m_r + m_r*other.m_i);
	}
	CxSqrt2Ring& operator*=(const CxSqrt2Ring& other) {
		CxSqrt2Ring s = (*this)*other;
		*this = std::move(s);
		return *this;
	}
	CxSqrt2Ring& operator/=(const CxSqrt2Ring& other) {
		CxSqrt2Ring s = (*this)/other;
		*this = std::move(s);
		return *this;
	}
	CxSqrt2Ring operator/(const CxSqrt2Ring& other) const {
		return (operator*(other.getConj())).operator/(Sqrt2Ring(other.re()*other.re() + other.im()*other.im()));
	}
	CxSqrt2Ring& operator*=(const Sqrt2Ring& i) {m_r *= i; m_i *= i; return *this;}
	CxSqrt2Ring& operator/=(const Sqrt2Ring& i) {m_r /= i; m_i /= i; return *this;}
	CxSqrt2Ring operator*(const Sqrt2Ring& i) const {CxSqrt2Ring ret(*this); return ret *= i;}
	CxSqrt2Ring operator/(const Sqrt2Ring& i) const {CxSqrt2Ring ret(*this); return ret /= i;}
	operator bool() const {return static_cast<bool>(m_r) || static_cast<bool>(m_i);}
};

#endif

#if USE_INT

using R = I;

#elif USE_SQRT2

using R = Sqrt2Ring;

#elif USE_CX

using R = CxSqrt2Ring;

#endif

class Matrix;

class ExteriorAlgebra {
	class BaseElement {
		std::vector<std::size_t> m_e;
		std::size_t m_hash = 0;
	
	public:
		BaseElement() =default;
		BaseElement(std::vector<std::size_t>&& e) noexcept;
		BaseElement& operator=(std::vector<std::size_t>&& e) noexcept;
		const std::vector<std::size_t>& getVector() const{return m_e;}
		std::size_t hash() const {return m_hash;}
		bool operator==(const BaseElement& other) const {return (m_hash==other.m_hash) && (m_e==other.m_e);}
		bool operator!=(const BaseElement& other) const {return !(*this == other);}
	};

	struct BaseElementHash {
		std::size_t operator()(const BaseElement& baseElement) const {return baseElement.hash();}
	};

	using MapType = std::unordered_map<BaseElement, R, BaseElementHash>;
	MapType m_e;
	static bool multElement(const std::pair<BaseElement, R>& p1, const std::pair<BaseElement, R>& p2,
			std::pair<BaseElement, R>& ret);
	ExteriorAlgebra& addElement(const std::pair<const BaseElement, R>& p);

	ExteriorAlgebra(MapType e) : m_e(e) {}
public:
	ExteriorAlgebra() =default;
	ExteriorAlgebra(const R& val);
	template<std::size_t N>	
	ExteriorAlgebra(std::initializer_list<std::pair<std::array<std::size_t, N>, R>> l);
	ExteriorAlgebra(const Matrix& m);

	std::size_t size() const {return m_e.size();}
	R getSingleElement() const {assert(m_e.size()<=1); if (m_e.size() == 1) return m_e.begin()->second; else return R(0);}

	ExteriorAlgebra& operator+=(const ExteriorAlgebra& other);
	ExteriorAlgebra& operator-=(const ExteriorAlgebra& other);
	ExteriorAlgebra operator+(const ExteriorAlgebra& other) const {ExteriorAlgebra ret(*this); return ret += other;}
	ExteriorAlgebra operator-(const ExteriorAlgebra& other) const {ExteriorAlgebra ret(*this); return ret -= other;}
	ExteriorAlgebra operator*(const ExteriorAlgebra& other) const;
	ExteriorAlgebra& operator*=(const ExteriorAlgebra& other) {return *this = (*this)*other;}
	ExteriorAlgebra& operator*=(const R& other);
	ExteriorAlgebra operator*(const R& other) const {ExteriorAlgebra ret(*this); return ret *= other;}
	bool operator==(const ExteriorAlgebra& other) const {return m_e==other.m_e;}
	bool operator!=(const ExteriorAlgebra& other) const {return m_e!=other.m_e;}
	operator bool() const {return size()!=0;}
};

class Matrix {
	std::vector<R> m_data;
	std::size_t m_size = 0;

	std::vector<R> svec() const;

public:
	Matrix() =default;
	explicit Matrix(std::size_t s) : m_data(s*s), m_size(s) {}
	Matrix(const Matrix& a, const Matrix& b, const Matrix& c, const Matrix& d);
	std::size_t size() const {return m_size;}
	R& operator()(std::size_t i, std::size_t j) {return m_data[i*m_size + j];}
	const R& operator()(std::size_t i, std::size_t j) const {return m_data[i*m_size + j];}
	Matrix operator*(const R& other) const;
	Matrix operator/(const R& other) const;
	Matrix operator*(const Matrix& other) const;
	Matrix operator+(const Matrix& other) const;
	Matrix operator-(const Matrix& other) const;
	Matrix operator-() const;
	Matrix multTranspose(const Matrix& other) const;
	Matrix kroneckerProduct(const Matrix& other) const;
	operator bool() const;
	R permProduct(const std::vector<std::size_t>& per) const;
	
	static Matrix getElementarySymMatrixKro(std::size_t i, std::size_t j, std::size_t s);
	static Matrix getElementarySymMatrixCurv(std::size_t i, std::size_t j, std::size_t s);
	static Matrix getIdentityMatrix(std::size_t s);
};

using A = ExteriorAlgebra;

class MatrixDet {
	std::vector<A> m_data;
	std::size_t m_size = 0;
	
public:
	std::size_t size() const {return m_size;}
	explicit MatrixDet(std::size_t s) : m_data(s*s), m_size(s) {}
	A& operator()(std::size_t i, std::size_t j) {return m_data[i*m_size + j];}
	const A& operator()(std::size_t i, std::size_t j) const {return m_data[i*m_size + j];}
	A determinant() const;
};

unsigned pairToIndex(std::pair<unsigned, unsigned> p);
std::pair<unsigned, unsigned> indexToPair(unsigned i);
unsigned orderedPairToIndex(std::pair<unsigned, unsigned> p, unsigned n);
std::pair<unsigned, unsigned> indexToOrderedPair(unsigned i, unsigned n);

template<std::size_t N>
ExteriorAlgebra::ExteriorAlgebra(std::initializer_list<std::pair<std::array<std::size_t, N>, R>> l) {
	for (const auto& p : l) {
		const R& val = p.second;
		if (val) {
			std::vector<std::size_t> v(p.first.cbegin(), p.first.cend());
			assert(std::is_sorted(v.begin(), v.end())); // Only support sorted generators
			for (std::size_t i=0; i<v.size()-1; ++i) { // Check for repeated generators
				if (v[i] == v[i+1]) continue;
			}
			auto i = m_e.insert(std::make_pair<ExteriorAlgebra::BaseElement, R>(std::move(v), R(val)));
			if (!i.second) {
				i.first->second += val;
			}
		}
	}
}


#if USE_SQRT2 or USE_CX

std::ostream& operator<<(std::ostream& out, const Sqrt2Ring val);

#endif

#if USE_CX

std::ostream& operator<<(std::ostream& out, const CxSqrt2Ring val);

#endif

std::ostream& operator<<(std::ostream& out, const Rational& r);

#endif

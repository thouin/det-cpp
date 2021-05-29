#include "det.hpp"

std::size_t pairToIndex(std::pair<std::size_t, std::size_t> p) {
	if (p.first > p.second) std::swap(p.first, p.second);
	return (p.second*(p.second+1))/2 + p.first;
}

std::pair<std::size_t, std::size_t> indexToPair(std::size_t i) {
	std::size_t x = 0;
	std::size_t y;
	while (((x+1)*(x+2))/2 <= i) ++x;
	y = i - x*(x+1)/2;
	return std::make_pair(y, x);
}

std::size_t orderedPairToIndex(std::pair<std::size_t, std::size_t> p, std::size_t n) {
	if (p.first == p.second)
		return p.first;
	else
		return (p.second*(p.second-1))/2 + p.first + n;
}

std::pair<std::size_t, std::size_t> indexToOrderedPair(std::size_t i, std::size_t n) {
	if (i < n)
		return std::make_pair(i, i);
	else {
		i -= n;
		std::size_t x = 0;
		std::size_t y;
		while ((x*(x+1))/2 <= i) ++x;
		y = i - x*(x-1)/2;
		return std::make_pair(y, x);
	}
}

#if USE_INT

namespace {
	ExteriorAlgebra getIJCurveForm(std::size_t i, std::size_t j, std::size_t s) {
		ExteriorAlgebra ret;
		for (std::size_t k=0; k<s; k++) {
			std::size_t r = 2*pairToIndex(std::make_pair(i, k));
			std::size_t s = 2*pairToIndex(std::make_pair(j, k))+1;
			if (r <= s) {
				std::array<std::size_t, 2> q({r, s});
				ret += ExteriorAlgebra({std::make_pair(q, R(1))});
			} else {
				std::array<std::size_t, 2> q({s, r});
				ret += ExteriorAlgebra({std::make_pair(q, R(-1))});
			}
		}
		return ret;
	}
}

#endif

ExteriorAlgebra getCurveForm(std::pair<std::size_t, std::size_t> p1, std::pair<std::size_t, std::size_t> p2, std::size_t s) {
#if USE_INT
	ExteriorAlgebra ret;
	
	std::size_t i, j;
	if (p1.first == p2.first) {
		i = p1.second;
		j = p2.second;
		ret += getIJCurveForm(i, j, s);
	}
	if (p1.first == p2.second) {
		i = p1.second;
		j = p2.first;
		ret += getIJCurveForm(i, j, s);
	}
	if (p1.second == p2.first) {
		i = p1.first;
		j = p2.second;
		ret += getIJCurveForm(i, j, s);
	}
	if (p1.second == p2.second) {
		i = p1.first;
		j = p2.first;
		ret += getIJCurveForm(i, j, s);
	}
	
	return ret;
#else
	Matrix a = (Matrix::getElementarySymMatrixCurv(p1.first, p1.second, s)*Matrix::getElementarySymMatrixCurv(p2.first, p2.second, s)
		- Matrix::getElementarySymMatrixCurv(p2.first, p2.second, s)*Matrix::getElementarySymMatrixCurv(p1.first, p1.second, s))
		.kroneckerProduct(Matrix::getIdentityMatrix(s));
	Matrix b = (Matrix::getElementarySymMatrixCurv(p1.first, p1.second, s)*Matrix::getElementarySymMatrixCurv(p2.first, p2.second, s))
		.kroneckerProduct(Matrix::getIdentityMatrix(s));
	Matrix c = -(Matrix::getElementarySymMatrixCurv(p2.first, p2.second, s)*Matrix::getElementarySymMatrixCurv(p1.first, p1.second, s))
		.kroneckerProduct(Matrix::getIdentityMatrix(s));
	Matrix d(s*(s+1)/2);
	
#if !USE_CX
	Matrix ret(d, b, c, d);
#else
	R i = R(0, 1);
	Matrix ret(d, b*i, c*i, d);
#endif
	return ExteriorAlgebra(ret);
#endif
}

int main() {
	std::size_t N = MAX_SIZE+1;
	std::cout << "Dimension number:" << std::endl;
	std::cin >> N;
	if (N>MAX_SIZE) {
		std::cout << "Invalid dimension number!" << std::endl;
		return 0;
	}
	const std::size_t M = (N*(N+1))/2;
	MatrixDet mat(M);
	std::pair<std::size_t, std::size_t> p1;
	std::pair<std::size_t, std::size_t> p2;
	
	#pragma omp parallel for collapse(2) private(p1, p2)
	for (std::size_t i=0; i<M; ++i) {
		for (std::size_t j=0; j<M; ++j) {
			p1 = indexToPair(i);
			p2 = indexToPair(j);
			mat(i, j) = getCurveForm(p1, p2, N);
		}
	}
	
	ExteriorAlgebra v = mat.determinant();
	
#if USE_INT
	//R res = ((N+3)/4)%2 ? -v.getSingleElement() : v.getSingleElement();
	R res = v.getSingleElement();
	std::cout << Rational(res, 1UL << (3*N*(N+1)/2)) << std::endl;
#elif USE_SQRT2
	R res = (N*(N+1)/2 + (N*(N+1)+1)/4)%2 ? -v.getSingleElement() : v.getSingleElement();
	std::cout << res / R(1UL << N*(N+1)/2) << std::endl;
#elif USE_CX
	R res = (N*(N+1)/2)%2 ? -v.getSingleElement() : v.getSingleElement();
	bool im = (N*(N+1)/2)%2;
	if (im) {
		if (res.re()) std::cout << "Complex determinant!" << std::endl;
		else {
			std::cout << res.im() / Sqrt2Ring(1UL << N*(N+1)/2) << std::endl;
		}
	} else {
		if (res.im()) std::cout << "Complex determinant!" << std::endl;
		else {
			std::cout << res.re() / Sqrt2Ring(1UL << N*(N+1)/2) << std::endl;
		}
	}
#endif
	return 0;
}

#include <mutex>

#include "det.hpp"

std::vector<R> Matrix::svec() const {
	std::vector<R> ret(size()*(size()+1)/2);
	std::size_t k=0;
	for (std::size_t i=0; i<size(); i++) {
		for (std::size_t j=i; j<size(); j++) {
			assert((*this)(i, j) == (*this)(j, i));
#if USE_INT
			ret[k++] = (*this)(i, j);
#else
			ret[k++] = (i==j) ? (*this)(i, j) : R(Sqrt2Ring(0, 1)) * (*this)(i, j);
#endif
		}
	}
	return ret;
}

Matrix::Matrix(const Matrix& a, const Matrix& b, const Matrix& c, const Matrix& d) : m_data(4*a.size()*a.size()), m_size(2*a.size()) {
	assert(a.size() == b.size() && b.size() == c.size() && c.size() == d.size());
	for (std::size_t i=0; i<a.size(); i++) {
		for (std::size_t j=0; j<a.size(); j++) {
			(*this)(i, j) = a(i, j);
			(*this)(i, j+a.size()) = b(i, j);
			(*this)(i+a.size(), j) = c(i, j);
			(*this)(i+a.size(), j+a.size()) = d(i, j);
		}
	}
}

Matrix Matrix::operator*(const R& other) const {
	Matrix ret(size());
	for (std::size_t i=0; i<size()*size(); i++) {
		ret.m_data[i] = m_data[i]*other;
	}
	return ret;
}

Matrix Matrix::operator/(const R& other) const {
	Matrix ret(size());
	for (std::size_t i=0; i<size()*size(); i++) {
		ret.m_data[i] = m_data[i]/other;
	}
	return ret;
}

Matrix Matrix::operator*(const Matrix& other) const {
	assert(size() == other.size());
	Matrix ret(size());
	for (std::size_t i=0; i<size(); i++) {
		for (std::size_t j=0; j<size(); j++) {
			for (std::size_t k=0; k<size(); k++) {
				ret(i, j) += (*this)(i, k)*other(k, j);
			}
		}
	}
	return ret;
}

Matrix Matrix::operator+(const Matrix& other) const {
	assert(size() == other.size());
	Matrix ret(size());
	for (std::size_t i=0; i<size()*size(); i++) {
		ret.m_data[i] = m_data[i] + other.m_data[i];
	}
	return ret;
}

Matrix Matrix::operator-(const Matrix& other) const {
	assert(size() == other.size());
	Matrix ret(size());
	for (std::size_t i=0; i<size()*size(); i++) {
		ret.m_data[i] = m_data[i] - other.m_data[i];
	}
	return ret;
}

Matrix Matrix::operator-() const {
	Matrix ret(size());
	for (std::size_t i=0; i<size()*size(); i++) {
		ret.m_data[i] = -m_data[i];
	}
	return ret;
}

Matrix Matrix::multTranspose(const Matrix& other) const {
	assert(size() == other.size());
	Matrix ret(size());
	for (std::size_t i=0; i<size(); i++) {
		for (std::size_t j=0; j<size(); j++) {
			for (std::size_t k=0; k<size(); k++) {
				ret(i, j) += (*this)(i, k)*other(j, k);
			}
		}
	}
	return ret;
}

Matrix Matrix::kroneckerProduct(const Matrix& other) const {
	assert(size() == other.size());
	Matrix ret(size()*(size()+1)/2);
	std::size_t k = 0;
	for (std::size_t i=0; i<size(); i++) {
		for (std::size_t j=i; j<size(); j++) {
			Matrix m(getElementarySymMatrixKro(i, j, size()));
			std::vector<R> v = ( (((*this)*m).multTranspose(other) + (other*m).multTranspose(*this))
#if !USE_INT
					 / R(2)
#endif
			).svec();
			for (auto it=v.begin(); it!=v.end(); it++) {
				ret.m_data[k++] = std::move(*it);
			}
		}
	}
	return ret;
}

Matrix::operator bool() const {
	for (std::size_t i=0; i<size()*size(); i++) {
		if (m_data[i]) return true;
	}
	return false;
}

Matrix Matrix::getElementarySymMatrixKro(std::size_t i, std::size_t j, std::size_t s) {
	Matrix ret(s);
	if (i==j) {
#if USE_INT
		ret(i, i) = R(2);
#else
		ret(i, i) = R(1);
#endif
	} else {
#if USE_INT
		ret(i, j) = R(1);
#else
		ret(i, j) = R(Sqrt2Ring(0, Rational(1, 2)));
#endif
		ret(j, i) = ret(i, j);
	}
	return ret;
}

Matrix Matrix::getElementarySymMatrixCurv(std::size_t i, std::size_t j, std::size_t s) {
	Matrix ret(s);
	if (i==j) {
		ret(i, i) = R(1);
	} else {
#if USE_INT
		ret(i, j) = R(1);
#else
		ret(i, j) = R(Sqrt2Ring(0, Rational(1, 2)));
#endif
		ret(j, i) = ret(i, j);
	}
	return ret;
}

Matrix Matrix::getIdentityMatrix(std::size_t s) {
	Matrix ret(s);
	for (std::size_t i=0; i<s; i++) {
		ret(i, i) = R(1);
	}
	return ret;
}

namespace {
	A doDet(const MatrixDet& mat, std::vector<std::size_t> j, int sign, A prevProd, std::mutex& outputMtx) {
#if PRINT_SKIP
		static int outputCounter = 0;
#endif
		
		if (mat.size() != j.size()) {
			A ret;
			A prod;
			std::array<std::size_t, MAX_MATRIX_SIZE> sortedJ;
			std::copy(j.begin(), j.end(), sortedJ.begin()); std::sort(sortedJ.begin(), sortedJ.begin()+j.size());
			std::size_t sortedIndex = 0;
			std::size_t k = 0;
			j.push_back(0);
			
			if (j.size() < PARALLEL_THRESHOLD) {
				std::vector<A> res(2*MAX_SIZE);
				std::size_t l = 0;
				
				for (std::size_t i=0; i<mat.size(); i++) {
					if (sortedIndex < j.size()-1 && i == sortedJ[sortedIndex]) {
						sortedIndex++;
						continue;
					}

					k++;
					if (mat(j.size()-1, i).size() == 0) continue; // If the element is zero we dont need to calculate
					prod = mat(j.size()-1, i)*prevProd;
					if (prod) {
						int newSign = (k % 2) ? sign : -sign;
						j[j.size()-1] = i;
						#pragma omp task shared(res, mat, outputMtx)
						res[l] = doDet(mat, j, newSign, prod, outputMtx);
						l++;
					}
				}
			
				#pragma omp taskwait
				for (std::size_t i=0; i<l; i++) ret+=res[i];
			} else {
				for (std::size_t i=0; i<mat.size(); i++) {
					if (sortedIndex < j.size()-1 && i == sortedJ[sortedIndex]) {
						sortedIndex++;
						continue;
					}

					k++;
					if (mat(j.size()-1, i).size() == 0) continue; // If the element is zero we dont need to calculate
					prod =  mat(j.size()-1, i)*prevProd;
					if (prod) {
						int newSign = (k % 2) ? sign : -sign;
						j[j.size()-1] = i;
						ret += doDet(mat, j, newSign, prod, outputMtx);
					}
				}
			}
			
			return ret;
		} else {
			A ret(A(R(sign))*prevProd);
			{
				std::lock_guard<std::mutex> l(outputMtx);
#if PRINT_SKIP
				outputCounter++;
				if ((outputCounter %= PRINT_SKIP+1) == 0) {
#endif
					for (std::size_t i = 0; i<mat.size(); i++)
						std::cout << j[i] << ' ';
					std::cout << ret.getSingleElement() << std::endl;
#if PRINT_SKIP
				}
#endif
			}
			return ret;
		}
	}
}

A MatrixDet::determinant() const {
	if (m_size==0) return A();
	std::mutex outputMtx;
	A ret;
	#pragma omp parallel
	{
		#pragma omp single nowait
		ret = doDet(*this, std::vector<std::size_t>(), 1, A(R(1)), outputMtx);
	}
	return ret;
}

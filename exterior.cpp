#include <cassert>
#include <utility>
#include <cstddef>
#include <climits>

#include "det.hpp"

namespace std {
	template<>
	struct hash<ExteriorAlgebra> {
		std::size_t operator()(const ExteriorAlgebra&) const {return 0;}
	};
}

ExteriorAlgebra::BaseElement::BaseElement(std::vector<std::size_t>&& e) noexcept : m_e(std::move(e)) {
	for (std::size_t i : e) {
		m_hash ^= 1 << (i % (sizeof(std::size_t)*CHAR_BIT));
	}
}

ExteriorAlgebra::BaseElement& ExteriorAlgebra::BaseElement::operator=(std::vector<std::size_t>&& e) noexcept {
	m_e = std::move(e);
	m_hash = 0;
	for (std::size_t i : e) {
		m_hash ^= 1 << (i % (sizeof(std::size_t)*CHAR_BIT));
	}
	return *this;
}

bool ExteriorAlgebra::multElement(const std::pair<typename ExteriorAlgebra::BaseElement, R>& p1,
		const std::pair<ExteriorAlgebra::BaseElement, R>& p2,
		std::pair<ExteriorAlgebra::BaseElement, R>& ret) {
	const std::vector<std::size_t>& vp1 = p1.first.getVector();
	const std::vector<std::size_t>& vp2 = p2.first.getVector();
	std::vector<std::size_t> vret;
	int sign = 1;
	std::size_t p1ToSkip = vp1.size();
	auto it1 = vp1.cbegin();
	auto it2 = vp2.cbegin();
	while (it1!=vp1.cend() && it2!=vp2.cend()) {
		if (*it1 == *it2)
			return false;
		else if (*it1 < *it2) {
			vret.push_back(*it1);
			--p1ToSkip;
			++it1;
		} else {
			vret.push_back(*it2);
			sign = (p1ToSkip%2==0) ? sign : -sign;
			++it2;
		}
	}
	for (;it1!=vp1.cend();++it1) vret.push_back(*it1);
	for (;it2!=vp2.cend();++it2) vret.push_back(*it2);
	ret.first = std::move(vret);
	ret.second = R(sign)*p1.second*p2.second;
	return true;
}

ExteriorAlgebra& ExteriorAlgebra::addElement(const std::pair<const ExteriorAlgebra::BaseElement, R>& p) {
	if (p.second ){
		auto i = m_e.insert(p);
		if (!i.second) { // If there is already the element
			i.first->second += p.second;
			if (!i.first->second) {
				m_e.erase(i.first);
			}
		}
	}
	return *this;
}

ExteriorAlgebra::ExteriorAlgebra(const R& val) {
	if (val)
		m_e.insert(std::make_pair(ExteriorAlgebra::BaseElement(), R(val)));
}

ExteriorAlgebra::ExteriorAlgebra(const Matrix& m) {
	for (std::size_t i=0; i<m.size(); i++) {
		for (std::size_t j=i; j<m.size(); j++) {
			assert(m(i, j) == -m(j, i));
			if (m(i, j)) {
				m_e.insert(std::make_pair<ExteriorAlgebra::BaseElement, R>(ExteriorAlgebra::BaseElement({i, j}), R(m(i, j))));
			}
		}
	}
}

ExteriorAlgebra& ExteriorAlgebra::operator+=(const ExteriorAlgebra& other) {
	for (const auto& p : other.m_e) {
		addElement(p);
	}
	return *this;
}

ExteriorAlgebra& ExteriorAlgebra::operator-=(const ExteriorAlgebra& other) {
	for (const auto& p : other.m_e) {
		auto p0(p); p0.second = -p0.second;
		addElement(p0);
	}
	return *this;
}

ExteriorAlgebra ExteriorAlgebra::operator*(const ExteriorAlgebra& other) const {
	ExteriorAlgebra ret;
	std::pair<ExteriorAlgebra::BaseElement, R> p;
	for (const auto& p0 : m_e) {
		for (const auto& p1 : other.m_e) {
			bool b = multElement(p0, p1, p);
			if (b)
				ret.addElement(p); // If the product is not zero
		}
	}
	return ret;
}

ExteriorAlgebra& ExteriorAlgebra::operator*=(const R& other) {
	for (auto& p : m_e) {
		p.second *= other;
	}
	return *this;
}

ExteriorAlgebra operator*(const R& a, const ExteriorAlgebra& b) {return b*a;}

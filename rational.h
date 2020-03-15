//
// Created by Михаил Лобанов on 27.10.2019.
//

#ifndef BIGINTEGER_RATIONAL_H
#define BIGINTEGER_RATIONAL_H

#include <iostream>
#include <vector>
#include <string>




class BigInteger {
private:

    bool _sign;

    std::vector<int> _abs_vector;

public:

    BigInteger() {
        _abs_vector = {0};
        _sign = true;
    }

    operator int() const {
        int ans = 0;
        int cur = 1;
        for (size_t i = 0; i < _abs_vector.size(); i++) {
            ans += cur * _abs_vector[i];
            cur *= 10;
        }
        if (!_sign)
            ans = -ans;
        return ans;
    }

    explicit operator bool() const {
        return !(_abs_vector.size() == 1 && _abs_vector[0] == 0);
    }

    BigInteger(int x) {
        BigInteger k;
        k = x;
        *this = k;
    }

    BigInteger(const std::string &s) {
        BigInteger ans;
        ans = s;
        _abs_vector = ans._abs_vector;
        _sign = ans._sign;
    }

    BigInteger &operator=(const std::string &other) {
        _abs_vector.clear();
        if (other.size() == 2 && other[0] == '-' && other[1] == '0') {
            _sign = true;
            _abs_vector = {0};
            return *this;
        }
        if (other[0] == '-')
            _sign = false;
        else
            _sign = true;
        if (!_sign) {
            _abs_vector.resize(other.size() - 1);
        } else {
            _abs_vector.resize(other.size());
        }
        int sign1;
        if (_sign)
            sign1 = 1;
        else
            sign1 = 0;
        for (int i = int(other.size()) - 1; i >= 1 - sign1; --i) {
            _abs_vector[other.size() - 1 - i] = int(other[i]) - int('0');
        }
        return *this;
    }

    BigInteger& operator=(const int &right1) {
        _abs_vector.clear();
        int right = right1;
        if (right < 0)
            _sign = false;
        else
            _sign = true;
        right = abs(right);
        while (right != 0) {
            _abs_vector.push_back(right % 10);
            right /= 10;
        }
        if (_abs_vector.size() == 0)
            _abs_vector.push_back(0);
        return *this;
    }

    std::string toString() const {
        std::string ans;
        ans.resize(_abs_vector.size());
        for (int i = int(_abs_vector.size()) - 1; i >= 0; i--) {
            ans[_abs_vector.size() - i - 1] = char(_abs_vector[i] + int('0'));
        }
        if (!_sign)ans = "-" + ans;
        return ans;
    }

    friend std::ostream &operator<<(std::ostream &stream, const BigInteger &a);

    friend std::istream &operator>>(std::istream &stream, BigInteger &a);

    bool operator==(const BigInteger &other) const {
        return _sign == other._sign && _abs_vector == other._abs_vector;
    }

    bool operator!=(const BigInteger &right) const {
        return !(*this == right);
    }

    bool operator<(const BigInteger &right) const {
        if (!_sign && right._sign)
            return true;
        if (_sign && !right._sign)
            return false;
        if (*this == right)
            return false;
        if (_abs_vector.size() < right._abs_vector.size()) {
            if (_sign)return true; else return false;
        }
        if (_abs_vector.size() > right._abs_vector.size()) {
            if (_sign)return false; else return true;
        }
        bool f = false;
        for (int i = int(_abs_vector.size()) - 1; i >= 0; --i) {
            if (_abs_vector[i] < right._abs_vector[i]) {
                f = true;
                break;
            }
            if (_abs_vector[i] > right._abs_vector[i]) {
                f = false;
                break;
            }
        }
        if (_sign)
            return f;
        else
            return !f;
    }

    bool operator<=(const BigInteger &right) const {
        return *this < right || *this == right;
    }

    bool operator>(const BigInteger &right) const {
        return right < *this;
    }

    bool operator>=(const BigInteger &right) const {
        return *this > right || *this == right;
    }

    BigInteger operator-() const {
        BigInteger ans = *this;
        if (ans == BigInteger(0))
            ans._sign = true;
        else
            ans._sign = !ans._sign;
        return ans;
    }

    BigInteger &operator++() {
        return *this += 1;
    }

    BigInteger operator++(int) {
        BigInteger ans = *this;
        *this += 1;
        return ans;
    }

    BigInteger &operator--() {
        return *this -= 1;
    }

    BigInteger operator--(int) {
        BigInteger ans = *this;
        *this -= 1;
        return ans;
    }

    BigInteger &operator+=(const BigInteger &right) {
        if (_sign == right._sign) {
            int it = 0;
            int transf = 0;
            while (size_t(it) < _abs_vector.size() || size_t(it) < right._abs_vector.size() || transf == 1) {
                if (size_t(it) >= _abs_vector.size()) {
                    _abs_vector.push_back(0);
                }
                int r = 0;
                if (size_t(it) < right._abs_vector.size())r = right._abs_vector[it];

                _abs_vector[it] = _abs_vector[it] + r + transf;
                transf = _abs_vector[it] / 10;
                _abs_vector[it] = _abs_vector[it] % 10;
                it++;
            }
            while (_abs_vector.size() > 1 && _abs_vector[_abs_vector.size() - 1] == 0)
                _abs_vector.pop_back();
        } else {
            bool min_el = false;
            if (abs1(right) < abs1(*this))min_el = true;
            int it = 0;
            int transf = 0;
            _sign = (!_sign && abs1(*this) <= abs1(right)) || (_sign && abs1(*this) >= abs1(right));
            while (size_t(it) < _abs_vector.size() || size_t(it) < right._abs_vector.size()) {
                if (size_t(it) >= _abs_vector.size()) {
                    _abs_vector.push_back(0);
                }
                int r = 0;
                if (size_t(it) < right._abs_vector.size())r = right._abs_vector[it];
                if (min_el)_abs_vector[it] = _abs_vector[it] - r + transf;
                else
                    _abs_vector[it] = r - _abs_vector[it] + transf;
                if (_abs_vector[it] < 0) {
                    _abs_vector[it] += 10;
                    transf = -1;
                } else transf = 0;
                it++;
            }
            while (_abs_vector.size() > 1 && _abs_vector[_abs_vector.size() - 1] == 0)
                _abs_vector.pop_back();
        }
        return *this;
    }

    BigInteger &operator-=(const BigInteger &right) {
        BigInteger n = -right;
        *this += -right;
        return *this;
    }

    BigInteger operator+(const BigInteger &right) const {
        BigInteger ans = *this;
        ans += right;
        return ans;
    }

    BigInteger operator-(const BigInteger &right) const {
        BigInteger ans = *this;
        ans -= right;
        return ans;
    }

    BigInteger &operator*=(const BigInteger &right) {
        if (*this == BigInteger(0) || right == BigInteger(0)) {
            return *this = BigInteger(0);
        }
        BigInteger ans = 0;
        for (size_t i = 0; i < right._abs_vector.size(); i++) {
            BigInteger p = 0;
            p._abs_vector.clear();
            for (size_t j = 0; j < i; j++)
                p._abs_vector.push_back(0);
            for (size_t j = 0; j < _abs_vector.size(); j++)
                p._abs_vector.push_back(_abs_vector[j]);
            for (int j = 0; j < right._abs_vector[i]; j++)
                ans += p;
        }
        if (right._sign == _sign)
            ans._sign = true;
        else
            ans._sign = false;
        *this = ans;
        return *this;
    }

private:
    BigInteger abs1(BigInteger a) {
        a._sign = true;
        return a;
    }

public:
    BigInteger absolute_value() {
        BigInteger ans = *this;
        ans._sign = true;
        return ans;
    }


    BigInteger operator*(const BigInteger &right) const {
        BigInteger ans = *this;
        ans *= right;
        return ans;
    }

    BigInteger &operator/=(const BigInteger &right) {
        if (*this == BigInteger(0)) {
            return *this;
        }
        BigInteger ans = BigInteger(0);
        ans._abs_vector.clear();
        BigInteger num = *this;
        BigInteger div = right;
        div._sign = true;
        num._sign = true;
        while (div <= num) {
            div *= 10;
        }
        while (div._abs_vector.size() >= right._abs_vector.size()) {
            int kol = 0;
            while (div <= num) {
                num -= div;
                kol++;
            }
            ans._abs_vector.push_back(kol);
            for (int i = 0; i < int(div._abs_vector.size()) - 1; i++) {
                std::swap(div._abs_vector[i], div._abs_vector[i + 1]);
            }
            div._abs_vector.pop_back();
        }
        for (size_t i = 0; i < ans._abs_vector.size() / 2; i++) {
            std::swap(ans._abs_vector[i], ans._abs_vector[ans._abs_vector.size() - 1 - i]);
        }
        while (ans._abs_vector.size() > 1 && ans._abs_vector[ans._abs_vector.size() - 1] == 0)
            ans._abs_vector.pop_back();
        if (right._sign == _sign)
            ans._sign = true;
        else
            ans._sign = false;
        if (ans._abs_vector.size() == 1 && ans._abs_vector[0] == 0)
            ans._sign = true;
        *this = ans;
        return *this;
    }

    BigInteger operator/(const BigInteger &right) const {
        BigInteger ans = *this;
        ans /= right;
        return ans;
    }

    BigInteger &operator%=(const BigInteger &right) {
        BigInteger k = *this / right;
        *this -= (k * right);
        return *this;
    }

    BigInteger operator%(const BigInteger &right) const {
        BigInteger ans = *this;
        ans %= right;
        return ans;
    }

    void ins(const unsigned int x) {
        _abs_vector.insert(_abs_vector.begin(), x, 0);
    }
};

inline BigInteger gcd(const BigInteger &a, const BigInteger &b) {
    if (a % b == BigInteger(0))
        return BigInteger(b);
    else {
        BigInteger aa = a;
        BigInteger bb = b;
        while (bb != BigInteger(0)) {
            aa %= bb;
            std::swap(aa, bb);
        }
        return aa;
    }

}

class Rational {
public:
    Rational() : _numerator(0), _denominator(1) {}

    Rational(const int a) : _numerator(a), _denominator(1) {}

    Rational(const BigInteger &a) : _numerator(a), _denominator(1) {}

    Rational(const int a, const int b) : _numerator(a), _denominator(b) {
        this->_make_irreducible();
    }

    Rational(const BigInteger &a, const BigInteger &b) : _numerator(a), _denominator(b) {
        this->_make_irreducible();
    }

    std::string toString() {
        std::string s;
        if (_denominator != BigInteger(1))
            s = _numerator.toString() + "/" + _denominator.toString();
        else
            s = _numerator.toString();
        return s;
    }

    Rational operator-() const {
        Rational ans = *this;
        ans._numerator = -ans._numerator;
        return ans;
    }

    Rational &operator+=(const Rational &other) {
        _numerator *= other._denominator;
        _numerator += _denominator * other._numerator;
        _denominator *= other._denominator;
        this->_make_irreducible();
        return *this;
    }

    Rational &operator-=(const Rational &other) {
        *this += -other;
        return *this;
    }

    Rational &operator*=(const Rational &other) {
        _numerator *= other._numerator;
        _denominator *= other._denominator;
        this->_make_irreducible();
        return *this;
    }

    Rational reverse() const {
        Rational ans;
        ans._numerator = _denominator;
        ans._denominator = _numerator;
        ans._make_irreducible();
        return ans;
    }

    Rational &operator/=(const Rational &other) {
        *this *= other.reverse();
        return *this;
    }

    Rational operator+(const Rational &other) const {
        Rational ans = *this;
        ans += other;
        return ans;
    }

    Rational operator-(const Rational &other) const {
        Rational ans = *this;
        ans -= other;
        return ans;
    }

    Rational operator*(const Rational &other) const {
        Rational ans = *this;
        ans *= other;
        return ans;
    }

    Rational operator/(const Rational &other) const {
        Rational ans = *this;
        ans /= other;
        return ans;
    }

    bool operator==(const Rational &other) const {
        return _numerator * other._denominator == _denominator * other._numerator;
    }

    bool operator<(const Rational &other) const {
        return _numerator * other._denominator < _denominator * other._numerator;
    }

    bool operator<=(const Rational &other) const {
        return *this < other || *this == other;
    }

    bool operator>(const Rational &other) const {
        return other < *this;
    }

    bool operator>=(const Rational &other) const {
        return other <= *this;
    }

    bool operator!=(const Rational &other) const {
        return !(*this == other);
    }

    std::string asDecimal(size_t precision = 0) const {
        if (*this < 0)
            return "-" + this->_abs_asDecimal(precision);
        else
            return this->_abs_asDecimal(precision);
    }

    explicit operator double() const {
        return std::stold(this->asDecimal(350));
    }


private:
    BigInteger _numerator;
    BigInteger _denominator;

    void _make_irreducible() {
        BigInteger reduce = gcd(_numerator, _denominator);
        _numerator /= reduce;
        _denominator /= reduce;
        if (_denominator < BigInteger(0)) {
            _denominator *= -1;
            _numerator *= -1;
        }
    }

    std::string _abs_asDecimal(size_t precision = 0) const {
        Rational tmp = *this;
        tmp._numerator = tmp._numerator.absolute_value();
        std::string ans = (tmp._numerator / tmp._denominator).toString();
        ans.push_back('.');
        tmp -= tmp._numerator / tmp._denominator;
        for (unsigned int i = 0; i < precision; ++i) {
            tmp *= Rational(10);
            ans += (tmp._numerator / tmp._denominator).toString();
            tmp -= tmp._numerator / tmp._denominator;
        }
        return ans;
    }
};

Rational operator/(const int a, const Rational &b) {
    return Rational(a) / b;
}



std::ostream &operator<<(std::ostream &stream, const BigInteger &a) {
    if (!a._sign)stream << '-';
    for (int i = int(a._abs_vector.size()) - 1; i >= 0; i--)
        stream << a._abs_vector[i];
    return stream;
}

std::istream &operator>>(std::istream &stream, BigInteger &a) {
    std::string s;
    stream >> s;
    a = s;
    return stream;
}

#endif //BIGINTEGER_RATIONAL_H

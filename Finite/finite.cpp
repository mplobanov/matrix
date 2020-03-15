//
// Created by Михаил Лобанов on 12.03.2020.
//
#include "../prime_check/is_prime.cpp"

template <unsigned int MOD>
class Finite {
private:
    unsigned long long n;
public:
    Finite() {
        n = 0;
    }

    Finite(int num) {
        while (num < 0)
            num += MOD;
        n = num % MOD;
    }

    Finite operator+(const Finite & other) {
        return (n + other.n) % MOD;
    }

    Finite operator-(const Finite & other) {
        return (max(n, other.n) - min(n, other.n)) % MOD;
    }

    Finite operator*(const Finite & other) {
        return (n * other.n) % MOD;
    }

    Finite degree(unsigned int deg) {
        if (deg == 0)
            return Finite(1);
        if (deg % 2 == 0) {
            Finite ans = degree(deg / 2);
            ans = ans * ans;
            return ans;
        }
        else {
            Finite ans = degree(deg - 1);
            ans.n = (ans.n * n) % MOD;
            return ans;
        }
    }

    unsigned int toInt() {
        return n;
    }

    Finite inv(){
        static_assert(is_prime<MOD>::value, "Num is not prime");
        return degree(MOD - 2);
    };

    Finite operator/(const Finite & other) {
        return this * other.inv();
    }

};

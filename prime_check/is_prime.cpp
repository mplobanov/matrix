//
// Created by Михаил Лобанов on 11.03.2020.
//
template <int n, int div>
struct is_prime_rec {
    static const bool value = (n % div != 0) && is_prime_rec<n, div - 1>::value;
};

template <int n>
struct is_prime_rec<n, 1> {
    static const bool value = true;
};

template <int n>
struct is_prime_rec<n, 0> {
    static const bool value = true;
};

template <int n>
struct is_prime {
    static const bool value = is_prime_rec<n, n - 1>::value;
};


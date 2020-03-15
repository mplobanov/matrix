//
// Created by Михаил Лобанов on 15.03.2020.
//
#include "rational.h"
#include <vector>

template <unsigned N, unsigned M, typename Field = Rational>
class Matrix {
private:
    std::vector<std::vector<Field>> array;
public:
    Matrix(int n, int m) {
        array = std::vector<std::vector<Field>>(n, std::vector<Field>(M, Field(0)));
    }

    Matrix(int n, int m, const std::vector<std::vector<Field>> & arr) : array(arr) {}

    Matrix operator+=(const Matrix<N, M, Field> & other) {
        for (int i = 0; i < array.size(); ++i) {
            for (int j = 0; j < array[0].size(); ++j) {
                array[i][j] += other.array[i][j];
            }
        }
        return this;
    }

    Matrix operator-() {
        Matrix ans = this;
        for (int i = 0; i < ans.size(); ++i) {
            for (int j = 0; j < ans[0].size(); ++j) {
                ans[i][j] = -ans.array[i][j];
            }
        }
        return ans;
    }

    Matrix operator-=(const Matrix<N, M, Field> & other){
        this += -other;
        return this;
    }

};


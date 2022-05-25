#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>
#include "windows.h"


class Operations {
public:
    template <typename T>
    static std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>>& mat) {
        std::vector<std::vector<T>> b(mat.size(), std::vector<T>(mat[0].size()));
        for (size_t i = 0; i < mat.size(); i++) {
            for (size_t j = 0; j < mat[0].size(); j++) {
                b[i][j] = mat[j][i];
            }
        }
        return b;
    }
    template <typename T>
    static std::vector<std::vector<T>> Multiply_Matrix(const std::vector<std::vector<T>>& lhs,
        const std::vector<std::vector<T>>& rhs) {
        size_t n = lhs.size(), m = lhs[0].size();
        std::vector<std::vector<T>> H(n, std::vector<T>(m));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                H[i][j] = 0;
                for (int t = 0; t < n; t++) {
                    H[i][j] += lhs[i][t] * rhs[t][j];
                }
            }
        }
        return H;
    }
};

std::vector<double> GaussMeth(std::vector<std::vector<double>>& matr, std::vector<double>& b) {
    double d, s;
    size_t n = matr.size();
    std::vector<double> res(n);
    for (size_t k = 0; k < n; ++k)
    {
        double max_el = DBL_MIN;
        size_t index = 0;
        for (size_t i = k; i < n; ++i)
        {
            if (max_el < std::fabs(matr[i][k]))
            {
                max_el = std::fabs(matr[i][k]);
                index = i;
            }
        }
        if (index != k)
        {
            for (size_t i = 0; i < n; ++i)
            {
                double temp = matr[index][i];
                matr[index][i] = matr[k][i];
                matr[k][i] = temp;
            }
            double t = b[index];
            b[index] = b[k];
            b[k] = t;
        }
        for (size_t j = k + 1; j < n; ++j)
        {

            d = matr[j][k] / matr[k][k];
            for (size_t i = k; i < n; ++i)
            {
                matr[j][i] = matr[j][i] - d * matr[k][i];
            }
            b[j] = b[j] - d * b[k];
        }
    }
    for (size_t k = n - 1; k >= 0; --k)
    {
        d = 0;
        for (size_t j = k; j < n; ++j)
        {
            s = matr[k][j] * (k == n - 1 ? 0 : res[j]);
            d += s;
        }

        res[k] = (b[k] - d) / matr[k][k];
        if (k == 0) {
            break;
        }
    }


    return res;
}

class Task {
private:
    double eps = 0.0001;
public:
    virtual void Execute() {}
    void SetEps(double value) {
        if (value < 0 || value >= 1) {
            return;
        }
        eps = value;
    }
    double GetEps() {
        return eps;
    }
};

class Newton : public Task {
private:
	double x0 = 1.25;
	double y0 = 0;
public:
    void Execute() override {
        std::cout << "Метод Ньютона" << std::endl;
        auto f1 = [](double x, double y) {
            return tan(x * y + 0.1) - pow(x, 2);
        };
        auto f2 = [](double x, double y) {
            return pow(x, 2) + 2.0 * pow(y, 2) - 1;
        };
        auto f1_derx = [](double x, double y) {
            return -2 * x + y / (pow(cos(0.1 + x * y), 2));
        };
        auto f1_dery = [](double x, double y) {
            return x / (pow(cos(0.1 + x * y), 2));
        };
        auto f2_derx = [](double x, double y) {
            return 2 * x;
        };
        auto f2_dery = [](double x, double y) {
            return 4 * y;
        };
        double x = x0, y = y0;
        std::vector<std::vector<double>> A(2, std::vector<double>(2));
        std::vector<double> F_(2);
        std::vector<double> sol = { DBL_MAX, DBL_MAX };
        while (max(fabs(sol.front()), fabs(sol.back())) > GetEps()) {
            A[0][0] = f1_derx(x, y);
            A[0][1] = f2_derx(x, y);
            A[1][0] = f1_dery(x, y);
            A[1][1] = f2_dery(x, y);
            F_[0] = f1(x, y);
            F_[1] = f2(x, y);
            sol = GaussMeth(A, F_);
            x -= sol.front();
            y -= sol.back();
        }
        std::cout << "Відповідь: " << std::endl;
        std::cout << "x = " << x << std::endl;
        std::cout << "y = " << y << std::endl;
	}
};
class Rotating_Jakobi : public Task {
private:
    std::vector<std::vector<double>> _matrix = {
        {2, 1, 0},
        {1, 2, 1},
        {0, 1, 2}
    };
    std::pair<size_t, size_t> findMax(const std::vector<std::vector<double>>& matrix) {
        double max = DBL_MIN;
        std::pair<size_t, size_t> indexes;
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[0].size(); ++j) {
                if (i == j) {
                    continue;
                }
                if (fabs(matrix[i][j]) > max) {
                    max = fabs(matrix[i][j]);
                    indexes = std::make_pair(i, j);
                }
            }
        }
        return indexes;
    }
    void fillU(std::vector<std::vector<double>>& U, const size_t& i, const size_t& j, const double& fi) {
        U[i][i] = cos(fi);
        U[i][j] = sin(fi);
        U[j][j] = cos(fi);
        U[j][i] = -sin(fi);
        for (size_t k = 0; k < U.size(); ++k) {
            if (k == i || k == j) {
                continue;
            }
            U[k][k] = 1.0;
        }
    }
    double GetSum(const std::vector<std::vector<double>>& A) {
        double sum = 0;
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = i + 1; j < A[0].size(); ++j) {
                sum += pow(A[i][j], 2);
            }
        }
        return sum*2;
    }
public:
    void Execute() override {
        int iteration = 1;
        std::cout << "Метод обертань(Якобі)" << std::endl;

        if (auto transoped = Operations::Transpose(_matrix); transoped != _matrix) {
            return;
        }
        std::vector<std::vector<double>> A = _matrix;
        std::vector<std::vector<std::vector<double>>> Uses;
        
        while (GetSum(A) >= GetEps()) {
            /*std::cout << "" << iteration << ".\nСума квадратів не діагональних елементів перетвореної матриці > " << GetSum(A) << std::endl;
            ++iteration;*/
            std::vector<std::vector<double>> U(_matrix.size(), std::vector<double>(_matrix[0].size()));
            double fi = 0;
            std::pair<size_t, size_t> indexes = findMax(A);
            double multiplier = 0;
            if (A[indexes.first][indexes.first] == A[indexes.second][indexes.second]) {
                multiplier = -M_PI / 2.0;
            }
            else {
                multiplier = atan(2.0 * A[indexes.first][indexes.second] / (A[indexes.first][indexes.first] - A[indexes.second][indexes.second]));
            }
            fi = 0.5 * multiplier;
            fillU(U, indexes.first, indexes.second, fi);
            Uses.push_back(U);

            auto U_trans = Operations::Transpose(U);

            auto temp = Operations::Multiply_Matrix(U_trans, A);
            A = Operations::Multiply_Matrix(temp, U);
        }
        std::cout << "Власні вектори: " << std::endl;
        std::vector<std::vector<double>> ans = Uses.front();
        for (size_t i = 1; i < Uses.size(); ++i) {
            ans = Operations::Multiply_Matrix(ans, Uses[i]);
        }
        for (size_t i = 0; i < ans.size(); ++i) {
            std::cout << "| ";
            for (size_t j = 0; j < ans[0].size(); ++j) {
                std::cout << ans[i][j] << " ";
            }
            std::cout << "|" << std::endl;
        }
        std::cout << "Власні числа: " << std::endl;
        for (size_t i = 0; i < A.size(); ++i) {
            std::cout << A[i][i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Сума квадратів не діагональних елементів матриці вихідної = " << GetSum(A) << std::endl;       
    }
};

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    std::cout << std::fixed << std::showpoint;
    //std::cout << std::setprecision(2);
    std::vector<std::shared_ptr<Task>> operations { std::shared_ptr<Task>(new Newton), std::shared_ptr<Task>(new Rotating_Jakobi) };

    for (auto& op : operations) {
        op->Execute();
    }

    return 0;
}

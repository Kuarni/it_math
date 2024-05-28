#include <iostream>
#include <string>
#include <cmath>

#define MAX_ERROR 0.0
#define ERROR_GRID_MULT 10

struct TaskData {
    TaskData(double lambda, size_t grid_size) : lambda(lambda), n(grid_size - 1),
                                                h(4 * M_PI / sqrt(lambda) / n),
                                                x(new double[grid_size]),
                                                y(new double[grid_size]) {
        if (h > sqrt(6 / lambda))
            std::cerr << "Warning: h > sqrt(6/lambda)" << std::endl;

        for (int i = 0; i < grid_size; i++)
            x[i] = h * i;
    }

    const double lambda;
    const size_t n;
    const double h;
    double *const x;
    double *const y;

    ~TaskData() {
        delete[] x;
        delete[] y;
    }
};

struct Tridiagonal {
    Tridiagonal(size_t n) : size(n - 1),
                            a(new double[size]()), b(new double[size]()),
                            c(new double[size]()), d(new double[size]()) {}

    const size_t size;
    double *const a;
    double *const b;
    double *const c;
    double *const d;

    void solve(double *y) {
        for (auto i = 2; i <= size; i++) {
            auto j = i - 1;
            auto w = a[j] / b[j - 1];
            b[j] -= w * c[j - 1];
            d[j] -= w * d[j - 1];
        }

        y[size] = d[size - 1] / b[size - 1];
        for (auto i = size - 1; i >= 1; i--) {
            auto j = i - 1;
            y[i] = (d[j] - c[j] * y[i + 1]) / b[j];
        }
    }

    ~Tridiagonal() {
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
    }
};

double metric_phi(TaskData &td, int i, int j) {
    if (i > j) std::swap(i, j);

    auto l = td.lambda, h = td.h;
    auto xb = td.x[i - 1], xm = td.x[i], xa = td.x[i + 1];

    if (i == j)
        return (l * pow(xm, 2) * xa - l * xm * pow(xa, 2) + (l * pow(xa, 3)) / 3 -
                l * pow(xm, 2) * xb + l * xm * pow(xb, 2) - (l * pow(xb, 3)) / 3 + xa - xb)
               / pow(h, 2);
    else if (i + 1 == j)
        return ((-1 / 6.0) * (-6 + l * pow((xm - xa), 2)) * (xm - xa)) / pow(h, 2);
    else return 0;
}

inline double metric_f_phi(TaskData &td, int i) {
    auto l = td.lambda, h = td.h;
    auto xb = td.x[i - 1], xm = td.x[i], xa = td.x[i + 1];

    return (2 * (-(xm - xa) * sqrt(l) * cos(xm * sqrt(l)) + sin(xm * sqrt(l)) -
                 sin(xa * sqrt(l))) + 2 * (-sqrt(l) * (xm - xb) * cos(xm * sqrt(l)) + sin(xm * sqrt(l)) -
                                           sin(sqrt(l) * xb))) / h;
}

void run_fem(TaskData &taskData) {
    Tridiagonal tridiagonal(taskData.n);

    for (int i = 1; i <= taskData.n - 1; i++) {
        auto j = i - 1;

        if (i - 1 >= 1)
            tridiagonal.a[j] = metric_phi(taskData, i - 1, i);

        if (i + 1 < taskData.n)
            tridiagonal.c[j] = metric_phi(taskData, i + 1, i);

        tridiagonal.b[j] = metric_phi(taskData, i, i);
        tridiagonal.d[j] = metric_f_phi(taskData, i);
    }

    tridiagonal.solve(taskData.y);
    taskData.y[0] = 0;
    taskData.y[taskData.n] = 0;
}

double phi(TaskData &td, size_t i, double x) {
    if (i == 0) {
        if (td.x[0] <= x && x <= td.x[1])
            return (td.x[1] - x) / td.h;
        else
            return 0;
    } else if (i == td.n) {
        if (td.x[td.n - 1] <= x && x <= td.x[td.n])
            return (x - td.x[td.n - 1]) / td.h;
        else
            return 0;
    } else {
        if ((td.x[i - 1] <= x) && (x <= td.x[i]))
            return (x - td.x[i - 1]) / td.h;
        else if ((td.x[i] <= x) && (x <= td.x[i + 1]))
            return (td.x[i + 1] - x) / td.h;
        else
            return 0;
    }
}

double calc(TaskData &td, double x) {
    size_t l = 0, r = td.n;

    while (r - l > 1) {
        auto mid = (l + r) / 2;
        if (x > td.x[mid])
            l = mid;
        else
            r = mid;
    }

    return td.y[l] * phi(td, l, x) + td.y[r] * phi(td, r, x);
}

inline double function(TaskData &td, double x) {
    return sin(sqrt(td.lambda) * x);
}

void check_error(TaskData &td) {
    double maxError = 0.0;
    auto newN = td.n * ERROR_GRID_MULT;
    auto newH = td.h / ERROR_GRID_MULT;
    for (size_t i = 0; i < newN; i++) {
        double x = i * newH;
        auto res = calc(td, x);
        auto real_res = function(td, x);
        double err = std::abs(res - real_res);

        maxError = std::max(maxError, err);
    }

    printf("max error: %0.8f\th^2: %0.8f\n", maxError, pow(td.h, 2));
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Using " << argv[0] << " lambda grid_size" << std::endl;
    }

    TaskData data(std::stod(argv[1]), std::stoul(argv[2]));

    run_fem(data);

    check_error(data);
}

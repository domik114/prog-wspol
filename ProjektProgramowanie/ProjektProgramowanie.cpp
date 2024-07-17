#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <thread>
#include <omp.h>

double function(double A, double B, double C, double x) {
    double ret = 0;
    for (int i = 0; i < 50; i++) {
        ret += (A * x * x * x * x * x * x * x * x * x * x) + (B * x * x * x * x * x * x * x * x * x) + (C * x * x * x * x * x * x * x * x);
    }
    return ret;
}

double function(double x) {
    //return (sin(A * x * x * x * x * x * x * x * x * x * x) + sin(B * x * x * x * x * x * x * x * x) + exp(C * x * x * x * x * x * x * x));
    return sin(x) * exp(-x);
}

double eulerMethod(double A, double B, double C, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double x = a + i * h;
        sum += function(A, B, C, x);
    }

    return h * sum;
}

double eulerMethodOPENMP(double A, double B, double C, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        sum += function(A, B, C, (a + i * h));
    }

    return h * sum;
}

void eulerMethodTHREAD(double A, double B, double C, double a, double b, int n, double& result) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double x = a + i * h;
        sum += function(A, B, C, x);
    }

    result = h * sum;
}

double gearMethod(double A, double B, double C, double a, double b, int n) {
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

    for (int i = 0; i < n; i++) {
        double x = a + i * h;

        double k1 = h * f[i];
        double k2 = h * function(A, B, C, x + h / 2.0);
        double k3 = h * function(A, B, C, x + h / 2.0);
        double k4 = h * function(A, B, C, x + h);

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        f[i + 1] = function(A, B, C, x + h);
    }

    double result = y[n];

    delete[] y;
    delete[] f;

    return result;
}

double gearMethodOPENMP(double A, double B, double C, double a, double b, int n) {
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

#pragma omp parallel for num_threads(16)
    for (int i = 0; i < n; i++) {
        double x = a + i * h;

        double k1 = h * f[i];
        double k2 = h * function(A, B, C, x + h / 2.0);
        double k3 = h * function(A, B, C, x + h / 2.0);
        double k4 = h * function(A, B, C, x + h);

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        f[i + 1] = function(A, B, C, x + h);   
    }

    double result = y[n];

    delete[] y;
    delete[] f;

    return result;
}

void gearMethodTHREAD(double A, double B, double C, double a, double b, int n, double& result2) {
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

    for (int i = 0; i < n; i++) {
        double x = a + i * h;

        double k1 = h * f[i];
        double k2 = h * function(A, B, C, x + h / 2.0);
        double k3 = h * function(A, B, C, x + h / 2.0);
        double k4 = h * function(A, B, C, x + h);

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        f[i + 1] = function(A, B, C, x + h);
    }

    double result = y[n];

    result2 = y[n];

    delete[] y;
    delete[] f;

}

double adamsBashforthMethod(double A, double B, double C, double a, double b, int n) {
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

    for (int i = 0; i < 3; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + h * f[i];
        f[i + 1] = function(A, B, C, x + h);
    }

    for (int i = 3; i < n; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + (h / 24.0) * (55 * f[i] - 59 * f[i - 1] + 37 * f[i - 2] - 9 * f[i - 3]);
        f[i + 1] = function(A, B, C, x + h);
    }

    double result = y[n];

    delete[] y;
    delete[] f;

    return result;
}

double adamsBashforthMethodOPENMP(double A, double B, double C, double a, double b, int n) {    
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

    for (int i = 0; i < 3; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + h * f[i];
        f[i + 1] = function(A, B, C, x + h);
    }

#pragma omp parallel for
    for (int i = 3; i < n; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + (h / 24.0) * (55 * f[i] - 59 * f[i - 1] + 37 * f[i - 2] - 9 * f[i - 3]);
        f[i + 1] = function(A, B, C, x + h);    
    }

    double result = y[n];

    delete[] y;
    delete[] f;

    return result;
}

void adamsBashforthMethodTHREAD(double A, double B, double C, double a, double b, int n, double &result) {
    double h = (b - a) / n;
    double* y = new double[n + 1];
    double* f = new double[n + 1];

    y[0] = function(A, B, C, a);
    f[0] = function(A, B, C, a);

    for (int i = 0; i < 3; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + h * f[i];
        f[i + 1] = function(A, B, C, x + h);
    }

    for (int i = 3; i < n; i++) {
        double x = a + i * h;

        y[i + 1] = y[i] + (h / 24.0) * (55 * f[i] - 59 * f[i - 1] + 37 * f[i - 2] - 9 * f[i - 3]);
        f[i + 1] = function(A, B, C, x + h);
    }

    result = y[n];

    delete[] y;
    delete[] f;
}

int main() {
    double a, b;
    int n;
    double A, B, C;
    A = 100000000;
    B = 100000000;
    C = 100000000;

    /*std::cout << "Podaj poczatek przedzialu: ";
    std::cin >> a;

    std::cout << "Podaj koniec przedzialu: ";
    std::cin >> b;

    std::cout << "Podaj liczbe podzialow: ";
    std::cin >> n;*/

    a = 1;
    b = 100000000;
    n = 10000000;
    /*
    {
        double start = omp_get_wtime();
        double result = eulerMethod(A, B, C, a, b, n);
        double end = omp_get_wtime();
        std::cout << "Wynik dla euler sekw: " << result << std::endl;
        std::cout << "Czas liczenia: " << end - start << " sekund\n" << std::endl;

        double start1 = omp_get_wtime();
        double result1 = eulerMethodOPENMP(A, B, C, a, b, n);
        double end1 = omp_get_wtime();
        std::cout << "Wynik dla euler OPENMP: " << result1 << std::endl;
        std::cout << "Czas liczenia: " << end1 - start1 << " sekund\n" << std::endl;

        double start3 = omp_get_wtime();
        const int num_threads3 = std::thread::hardware_concurrency();
        std::vector<std::thread> threads3(num_threads3);
        std::vector<double> results3(num_threads3);
        int chunk_size3 = n / num_threads3;
        int remainder3 = n % num_threads3;
        int index3 = 0;

        for (int i = 0; i < num_threads3; i++) {
            int chunk3 = chunk_size3 + (i < remainder3 ? 1 : 0);
            threads3[i] = std::thread(eulerMethodTHREAD, A, B, C, a + index3, a + index3 + chunk3, chunk3, std::ref(results3[i]));
            index3 += chunk3;
        }

        for (int i = 0; i < num_threads3; i++) {
            threads3[i].join();
        }

        double final_result3 = 0.0;
        for (int i = 0; i < num_threads3; i++) {
            final_result3 += results3[i];
        }

        double end3 = omp_get_wtime();

        std::cout << "Wynik dla euler THREAD: " << final_result3 << std::endl;
        std::cout << "Czas liczenia: " << end3 - start3 << " sekund\n\n" << std::endl;
    }*/
    /*
    {
        double start = omp_get_wtime();
        double result = gearMethod(A, B, C, a, b, n);
        double end = omp_get_wtime();
        std::cout << "Wynik dla gear sekwencyjnie: " << result << std::endl;
        std::cout << "Czas liczenia: " << end - start << " sekund\n" << std::endl;

        double start1 = omp_get_wtime();
        double result1 = gearMethodOPENMP(A, B, C, a, b, n);
        double end1 = omp_get_wtime();
        std::cout << "Wynik dla gear OPENMP: " << result1 << std::endl;
        std::cout << "Czas liczenia: " << end1 - start1 << " sekund\n" << std::endl;

        double start3 = omp_get_wtime();
        const int num_threads3 = std::thread::hardware_concurrency();
        std::vector<std::thread> threads3(num_threads3);
        std::vector<double> results3(num_threads3);
        int chunk_size3 = n / num_threads3;
        int remainder3 = n % num_threads3;
        int index3 = 0;

        for (int i = 0; i < num_threads3; i++) {
            int chunk3 = chunk_size3 + (i < remainder3 ? 1 : 0);
            threads3[i] = std::thread(gearMethodTHREAD, A, B, C, a + index3, a + index3 + chunk3, chunk3, std::ref(results3[i]));
            index3 += chunk3;
        }

        for (int i = 0; i < num_threads3; i++) {
            threads3[i].join();
        }

        double final_result3 = 0.0;
        for (int i = 0; i < num_threads3; i++) {
            final_result3 += results3[i];
        }

        double end3 = omp_get_wtime();

        std::cout << "Wynik dla gear THREAD: " << final_result3 << std::endl;
        std::cout << "Czas liczenia: " << end3 - start3 << " sekund\n\n" << std::endl;
    }*/
    
    {
        double start = omp_get_wtime();
        double result = adamsBashforthMethod(A, B, C, a, b, n);
        double end = omp_get_wtime();
        std::cout << "Wynik dla adams-bashforth sekwencyjnie: " << result << std::endl;
        std::cout << "Czas liczenia: " << end - start << " sekund\n" << std::endl;

        double start2 = omp_get_wtime();
        double result2 = adamsBashforthMethodOPENMP(A, B, C, a, b, n);
        double end2 = omp_get_wtime();
        std::cout << "Wynik dla adams-bashforth OPENMP: " << result2 << std::endl;
        std::cout << "Czas liczenia: " << end2 - start2 << " sekund\n" << std::endl;

        double start3 = omp_get_wtime();
        const int num_threads3 = std::thread::hardware_concurrency();
        std::vector<std::thread> threads3(num_threads3);
        std::vector<double> results3(num_threads3);
        int chunk_size3 = n / num_threads3;
        int remainder3 = n % num_threads3;
        int index3 = 0;

        for (int i = 0; i < num_threads3; i++) {
            int chunk3 = chunk_size3 + (i < remainder3 ? 1 : 0);
            threads3[i] = std::thread(adamsBashforthMethodTHREAD, A, B, C, a + index3, a + index3 + chunk3, chunk3, std::ref(results3[i]));
            index3 += chunk3;
        }

        for (int i = 0; i < num_threads3; i++) {
            threads3[i].join();
        }

        double final_result3 = 0.0;
        for (int i = 0; i < num_threads3; i++) {
            final_result3 += results3[i];
        }

        double end3 = omp_get_wtime();

        std::cout << "Wynik dla adams-bashforth THREAD: " << final_result3 << std::endl;
        std::cout << "Czas liczenia: " << end3 - start3 << " sekund\n\n" << std::endl;
    }
    
    return 0;
}

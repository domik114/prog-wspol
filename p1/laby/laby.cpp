#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdio.h>

double f(double x, double A, double B, double C, double D) {
    //return(A * (x * x * x) + B * (x * x) + C * x + D);
    return sin(x) * exp(-x);
}

void calka(double A, double B, double C, double D, double x1, double x2, int n) {
    double h = (x2 - x1) / n;
    double s = 0;

    for (int i = 1; i < n; i++) {
        s += f((x1 + i * h), A, B, C, D);
    }

    s = (s + (f(x1, A, B, C, D) + f(x2, A, B, C, D)) / 2) * h;

    std::cout << "Wartosc calki wynosi: " << s << std::endl;
}

void calka_rown(double A, double B, double C, double D, double x1, double x2, int n) {  
    double h = (x2 - x1) / n;
    double s = 0;

#pragma omp parallel for reduction(+:s)
    for (int i = 1; i < n; i++) {
        s += f((x1 + i * h), A, B, C, D);
    }

    s = (s + (f(x1, A, B, C, D) + f(x2, A, B, C, D)) / 2) * h;

    std::cout << "Wartosc calki wynosi: " << s << std::endl;
}

int main()
{
    double A, B, C, D, x1, x2, h;
    int n;
    double s = 0;

    /*std::cout << "Podaj A: "; std::cin >> A;
    std::cout << "Podaj B: "; std::cin >> B;
    std::cout << "Podaj C: "; std::cin >> C;
    std::cout << "Podaj D: "; std::cin >> D;

    std::cout << "Podaj poczatek: "; std::cin >> x1;
    std::cout << "Podaj koniec: "; std::cin >> x2;
    std::cout << "Podaj n: "; std::cin >> n;*/

    A = 10;
    B = 10;
    C = 10;
    D = 10;

    x1 = 1;
    x2 = 100000000;
    n = 1000000000;

    double start = omp_get_wtime();
    calka(A, B, C, D, x1, x2, n);
    double stop = omp_get_wtime();

    double start1 = omp_get_wtime();
    calka_rown(A, B, C, D, x1, x2, n);
    double stop1 = omp_get_wtime();

    std::cout << "\nCzas sekwencyjny: " << stop - start << std::endl;
    std::cout << "Czas rownolegly: " << stop1 - start1 << std::endl;

    return 0;
}
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <string>
#include <iomanip>

using namespace std;

vector<vector<double>> dodawanie(const vector<vector<double>>& A, const vector<vector<double>>& B, double size) {
    vector<vector<double>> matrix(size, vector<double>(size));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = A[i][j] + B[i][j];
        }
    }

    return matrix;
}

vector<vector<double>> mnozenie(const vector<vector<double>>& A, const vector<vector<double>>& B, double size) {
    vector<vector<double>> matrix(size, vector<double>(size, 0));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                matrix[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return matrix;
}

vector<vector<double>> transpozycja(const vector<vector<double>>& A, double size) {
    vector<vector<double>> matrix(size, vector<double>(size));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[j][i] = A[i][j];
        }
    }

    return matrix;
}
vector<vector<double>> dodawanie_rownolegle(const vector<vector<double>>& A, const vector<vector<double>>& B, int size) {
    vector<vector<double>> matrix(size, vector<double>(size));

#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = A[i][j] + B[i][j];
        }
    }

    return matrix;
}

vector<vector<double>> mnozenie_rownolegle(const vector<vector<double>>& A, const vector<vector<double>>& B, int size) {
    vector<vector<double>> matrix(size, vector<double>(size, 0));

#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                matrix[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return matrix;
}

vector<vector<double>> transpozycja_rownolegle(const vector<vector<double>>& A, int size) {
    vector<vector<double>> matrix(size, vector<double>(size));

#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[j][i] = A[i][j];
        }
    }

    return matrix;
}

int main()
{
    ifstream fileA("A.csv");
    vector<vector<double>> A;

    string line;
    double size = 0;

    string sline;
    getline(fileA, sline);
    size = stoi(sline);

    while (getline(fileA, line)) {
        vector<double> row;
        size_t pos = 0;
        string t;

        while ((pos = line.find(";")) != string::npos) {
            t = line.substr(0, pos);
            row.push_back(stod(t));
            line.erase(0, pos + 1);
        }

        A.push_back(row);
    }

    ifstream fileB("B.csv");
    vector<vector<double>> B;

    double sizeB = 0;
    getline(fileB, sline);
    sizeB = stoi(sline);

    while (getline(fileB, line)) {
        vector<double> row;
        size_t pos = 0;
        string t;

        while ((pos = line.find(";")) != string::npos) {
            t = line.substr(0, pos);
            row.push_back(stod(t));
            line.erase(0, pos + 1);
        }

        B.push_back(row);
    }

    double start = omp_get_wtime();
    vector<vector<double>> licz = mnozenie(transpozycja(A, size), dodawanie(A, B, size), size);
    vector<vector<double>> C = mnozenie(licz, B, size);
    double end = omp_get_wtime();
    printf("Czas wykonania sekwencynie: %f sekund\n", end - start);

    double start_r = omp_get_wtime();
    vector<vector<double>> licz2 = mnozenie_rownolegle(transpozycja_rownolegle(A, size), dodawanie_rownolegle(A, B, size), size);
    vector<vector<double>> C2 = mnozenie_rownolegle(licz2, B, size);
    double end_r = omp_get_wtime();
    printf("Czas wykonania rownolegle: %f sekund\n", end_r - start_r);

    char name[100];
    sprintf_s(name, "C_%0.4f_%0.4f.csv", end - start, end_r - start_r);
    ofstream fileC(name);
    fileC << size << endl;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j < size - 1)
                fileC << fixed << setprecision(4) << C[i][j] << ";";
            else
                fileC << fixed << setprecision(4) << C[i][j];
        }

        fileC << endl;
    }

    fileA.close();
    fileB.close();
    fileC.close();

    return 0;
}
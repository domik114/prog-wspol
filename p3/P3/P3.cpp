#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <string>
#include <omp.h>

using namespace std;

vector<vector<double>> load(string name) {
    ifstream plik(name);
    vector<vector<double>> matrix;

    string sline;
    getline(plik, sline);
    int ilosc_rownan = stoi(sline);

    while (getline(plik, sline)) {
        vector<double> wiersz;
        stringstream ss(sline);
        string t;

        while (getline(ss, t, ';')) {
            wiersz.push_back(stod(t));
        }

        matrix.push_back(wiersz);
    }

    return matrix;
}

void save(vector<double>& wektor, double sekw, double rown) {
    char file_name[100];
    sprintf_s(file_name, "C_%.6f_%.6f.csv", sekw, rown);

    ofstream plik(file_name);
    plik << wektor.size() << endl;

    for (int i = 0; i < wektor.size(); ++i) {
        plik << fixed << setprecision(4) << wektor[i];
        if (i != wektor.size() - 1) {
            plik << ";";
        }
    }
    plik << endl;
}

void gauss_rown(vector<vector<double>>& matrix, vector<double>& wektor_x) {

    int n = matrix.size();

    for (int r = 0; r < n - 1; r++) {
#pragma omp parallel for
        for (int i = r + 1; i < n; i++) {
            for (int j = r + 1; j < n + 1; j++) 
                matrix[i][j] -= (matrix[i][r] / matrix[r][r]) * matrix[r][j];
            
            //matrix[i][r] = 0;
        }
    }

    int s = n - 1;

    wektor_x[s] = matrix[s][n] / matrix[s][s];

    for (int i = n - 1; i >= 0; i--) {
        double s = 0;

#pragma omp parallel for reduction(+:s)
        for (int r = i + 1; r < n; r++) 
            s += matrix[i][r] * wektor_x[r];
        
        wektor_x[i] = (matrix[i][n] - s) / matrix[i][i];
    }

}

void gauss(vector<vector<double>>& matrix, vector<double>& wektor_x) {

    int n = matrix.size();

    for (int r = 0; r < n - 1; r++) {
        for (int i = r + 1; i < n; i++) {
            for (int j = r + 1; j < n + 1; j++) 
                matrix[i][j] -= (matrix[i][r] / matrix[r][r]) * matrix[r][j];
            
            //matrix[i][r] = 0;
        }
    }

    int s = n - 1;

    wektor_x[s] = matrix[s][n] / matrix[s][s];

    for (int i = n - 1; i >= 0; i--) {
        double s = 0;
        for (int r = i + 1; r < n; r++) 
            s += matrix[i][r] * wektor_x[r];
        
        wektor_x[i] = (matrix[i][n] - s) / matrix[i][i];
    }
}


int main()
{
    double start, end, start_rownolegle, end_rownolegle;
    vector<vector<double>> matrix = load("C500.csv");
    int n = matrix.size();
    vector<double> wynik(n);
    vector<double> wynik2(n);

    start = omp_get_wtime();
    //gauss(matrix, wynik);
    end = omp_get_wtime();
    printf("Czas wykonania sekwencynie: %f sekund\n", end - start);

    start_rownolegle = omp_get_wtime();
    gauss_rown(matrix, wynik2);
    end_rownolegle = omp_get_wtime();
    printf("czas wykonania rownolegle: %f sekund\n", end_rownolegle - start_rownolegle);

    double czas_sekw = end - start;
    double czas_rown = end_rownolegle - start_rownolegle;

    save(wynik2, czas_sekw, czas_rown);

    return 0;
}
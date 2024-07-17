﻿#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <ctime>

using namespace std;

vector<vector<double>> mnozenie(vector<vector<double>>& A, vector<vector<double>>& B, int n, int m, int p) {
    vector<vector<double>> C(n, vector<double>(p, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

int main() {
    int rank, size;

    int n = 500; int m = 500; int p = 500;

    double czas_sekw;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<vector<double>> A, B, C;

    if (rank == 0) {        
        ifstream file_A("C500.csv");
        A = vector<vector<double>>(n, vector<double>(m));

        for (int i = 0; i < n; i++) {
            string line;
            getline(file_A, line);
            stringstream ss(line);
            for (int j = 0; j < m; j++) {
                string value;
                getline(ss, value, ';');
                A[i][j] = stod(value);
            }
        }

        file_A.close();    

        ifstream file_B("c500.csv");
        B = vector<vector<double>>(m, vector<double>(p));

        for (int i = 0; i < m; i++) {
            string line;
            getline(file_B, line);
            stringstream ss(line);
            for (int j = 0; j < p; j++) {
                string value;
                getline(ss, value, ';');
                B[i][j] = stod(value);
            }
        }

        file_B.close();        

        double start_sekw = MPI_Wtime();
        vector<vector<double>> C = mnozenie(A, B, n, m, p);

        double end_sekw = MPI_Wtime();
        czas_sekw = end_sekw - start_sekw;
    }

    double start_mpi = MPI_Wtime();

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);


    if (rank != 0) {
        A = vector<vector<double>>(n, vector<double>(m));
        B = vector<vector<double>>(m, vector<double>(p));
    }

    C = vector<vector<double>>(n, vector<double>(p));

    for (int i = 0; i < n; i++) 
        MPI_Bcast(A[i].data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < m; i++) 
        MPI_Bcast(B[i].data(), p, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

    int process_rows = n / size;
    int start_row = rank * process_rows;
    int end_row = (rank == size - 1) ? n : start_row + process_rows;

    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    if (rank == 0) {
        for (int k = 1; k < size; k++) {
            start_row = k * process_rows;
            end_row = (k == size - 1) ? n : start_row + process_rows;
            for (int i = start_row; i < end_row; i++) 
                MPI_Recv(C[i].data(), p, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);            
        }
    }
    else
    {
        for (int i = start_row; i < end_row; i++) 
            MPI_Send(C[i].data(), p, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);        
    }

    double end_mpi = MPI_Wtime();
    double czas_mpi = end_mpi - start_mpi;

    if (rank == 0) {
        char name[100];
        sprintf_s(name, "C_%0.4f_%0.4f.csv", czas_sekw, czas_mpi);
        ofstream plik_c(name);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                if (j < p - 1)
                    plik_c << fixed << setprecision(4) << C[i][j] << ";";
                else
                    plik_c << fixed << setprecision(4) << C[i][j];
            }
            plik_c << endl;
        }
        plik_c.close();
    }

    MPI_Finalize();

    return 0;
}


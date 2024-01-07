#include<iostream>
#include <mpi.h>
using namespace std;


int main(int argc, char** argv) {
   
    int size = 1000;
    int comm_rank, comm_size, thread, rows;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Status status;


    rows = size / comm_size;
    double* matr1 = new double[size * size];
    double* matr2 = new double[size * size];
    double* result = new double[size * size];

    if (comm_rank == 0) {

       for (int i = 0; i < size * size; ++i) {
            matr1[i] = i;
            matr2[i] = i;
            result[i] = 0;
        }
        
        start_time = MPI_Wtime();

        for (thread = 1; thread < comm_size; thread++) {

            MPI_Send(&matr1[thread * (rows * size) ], rows * size, MPI_DOUBLE, thread, 1, MPI_COMM_WORLD);
            MPI_Send(&matr2[0], size * size, MPI_DOUBLE, thread, 2, MPI_COMM_WORLD);
        }
    }

    else {

        MPI_Recv(&matr1[comm_rank * (rows * size) ], rows  * size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
       MPI_Recv(&matr2[0], size * size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    }


    for (int i = comm_rank * rows; i < comm_rank * rows+rows; i++) {  
        for (int j = 0; j < size; j++) {
            result[i * size + j] = 0;
            for (int t = 0; t < size; t++) {
               result[i * size + j] += matr1[i * size + t] * matr2[t * size + j];
            }
        }
    }


    if (comm_rank != 0) {
        MPI_Send(&result[comm_rank * (rows * size) ], rows  * size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    }
    else {

        for (thread = 1; thread < comm_size; thread++) {
            MPI_Recv(&result[thread * (rows*size) ], rows * size, MPI_DOUBLE, thread
              , 4, MPI_COMM_WORLD, &status);
        }

        end_time = MPI_Wtime();
        cout << end_time-start_time << endl;
    }
        

    MPI_Finalize();

    
    delete[] matr1;
    delete[] matr2;
    delete[] result;
    return 0;
}

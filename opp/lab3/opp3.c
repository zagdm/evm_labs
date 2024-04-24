#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>

#define Dim_count 2
#define Dim_X 0
#define Dim_Y 1
#define Grid_Size_max 16


void matrix_mult(double* A, double* B, double* C, int n1, int n2, int n3) {
    for (int i = 0; i < n1; ++i) {
        double* c = C + i * n3;
        for (int k = 0; k < n3; ++k) {
            c[k] = 0;
        }
        for (int j = 0; j < n2; ++j) {
            const double* b = B + j * n3;
            double a = A[i * n2 + j];
            for (int k = 0; k < n3; ++k) {
                c[k] += a * b[k];
            }
        }
    }
}

    
int correct_matrix_mult(double* A, double* B, double* C, double* tmp, int n1, int n2, int n3) {
    for (int i = 0; i < n1; ++i) {
        for (int k = 0; k < n3; k++) {
            tmp[k] = 0;
        }

        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                tmp[k] += A[i * n2 + j] * B[j * n3 + k];
            }
        }

        for (int k = 0; k < n3; k++) {
            if (tmp[k] != C[i * n3 + k]) return 0;
        }
    }
    return 1;
}

int main_exit(int exit_code, char* s_exit) {
    fputs(s_exit, stderr);
    MPI_Finalize();
    return exit_code;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        return main_exit(1, "Incorrect arguments\n");
    }
    const size_t n1 = strtoll(argv[1], NULL, 10);
    const size_t n2 = strtoll(argv[2], NULL, 10);
    const size_t n3 = strtoll(argv[3], NULL, 10);
    if (errno) {
        return main_exit(1, "Incorrect arguments\n");
    }


    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size > Grid_Size_max) {
        return main_exit(2, "Proc count incorrect\n");
    }

    int dims[Dim_count];
    dims[Dim_X] = dims[Dim_Y] = 0;
    MPI_Dims_create(size, Dim_count, dims);
    if (size != dims[Dim_X] * dims[Dim_Y]) {
        return main_exit(3, "Proc count incorrect\n");
    }
    if ((n1 % dims[Dim_X]) || (n3 % dims[Dim_Y])) {
        return main_exit(4, "Data size incorrect\n");
    }

    int periods[Dim_count];
    periods[Dim_X] = periods[Dim_Y] = 0;
    MPI_Comm gridComm;
    MPI_Cart_create(MPI_COMM_WORLD, Dim_count, dims, periods, 1, &gridComm);
    
    int coords[Dim_count];
    MPI_Comm columnComm;
    MPI_Comm rowComm;
    MPI_Comm_rank(gridComm, &rank);
    MPI_Cart_coords(gridComm, rank, Dim_count, coords);
    MPI_Comm_split(gridComm, coords[Dim_Y], coords[Dim_X], &columnComm);
    MPI_Comm_split(gridComm, coords[Dim_X], coords[Dim_Y], &rowComm);


    double* base_mem = (rank != 0) ? NULL : (double*)malloc((n1 * n2 + n3 * n2 + n1 * n3) * sizeof(double));
    double* base_A = (rank != 0) ? NULL : base_mem;
    double* base_B = (rank != 0) ? NULL : base_mem + n1 * n2;
    double* base_C = (rank != 0) ? NULL : base_mem + (n1 + n3) * n2;
    if (rank == 0) {
        if (base_mem == NULL) {
            return main_exit(5, "Out of memory");
        }
        srand(n1 * n2 * n3);
        for (size_t i = 0; i < (n1 + n3) * n2; i++) {
            base_mem[i] = rand();
        }
    }


    double start = MPI_Wtime();
    size_t row_size = n1 / dims[Dim_X];
    size_t column_size = n3 / dims[Dim_Y];
    int mpi_double_size;
    MPI_Type_size(MPI_DOUBLE, &mpi_double_size);

    double* part_mem = (double*)malloc((row_size * n2 + column_size * n2 + row_size * column_size) * sizeof(double));
    if (part_mem == NULL) {
        return main_exit(5, "Out of memory");
    }
    double* part_A = part_mem;
    double* part_B = part_mem + row_size * n2;
    double* part_C = part_mem + (row_size + column_size) * n2;


    if (coords[Dim_Y] == 0) {
        MPI_Scatter(base_A, row_size * n2, MPI_DOUBLE, part_A, row_size * n2, MPI_DOUBLE, 0, columnComm);
    }
    MPI_Bcast(part_A, row_size * n2, MPI_DOUBLE, 0, rowComm);


    if (coords[Dim_X] == 0) {
        MPI_Datatype column_type;
        MPI_Datatype compress_column_type;
        MPI_Type_vector(n2, column_size, n3, MPI_DOUBLE, &column_type);
        MPI_Type_create_resized(column_type, 0, column_size * mpi_double_size, &compress_column_type);
        MPI_Type_commit(&compress_column_type);
        MPI_Scatter(base_B, 1, compress_column_type, part_B, column_size * n2, MPI_DOUBLE, 0, rowComm);
        
        MPI_Type_free(&column_type);
        MPI_Type_free(&compress_column_type);
    }
    MPI_Bcast(part_B, column_size * n2, MPI_DOUBLE, 0, columnComm);


    matrix_mult(part_A, part_B, part_C, row_size, n2, column_size);


    MPI_Datatype matrix_rec_type;
    MPI_Datatype compress_matrix_rec_type;
    MPI_Type_vector(row_size, column_size, n3, MPI_DOUBLE, &matrix_rec_type);
    MPI_Type_create_resized(matrix_rec_type, 0, column_size * mpi_double_size, &compress_matrix_rec_type);
    MPI_Type_commit(&compress_matrix_rec_type);

    int size_table[Grid_Size_max];
    int shift_table[Grid_Size_max];
    for (int i = 0; i < size; i++) {
        size_table[i] = 1;
        shift_table[i] = i % dims[Dim_Y] + (i / dims[Dim_Y]) * dims[Dim_Y] * row_size;
    }

    MPI_Gatherv(part_C, column_size* row_size, MPI_DOUBLE, base_C, size_table, shift_table, compress_matrix_rec_type, 0, MPI_COMM_WORLD);
    MPI_Type_free(&matrix_rec_type);
    MPI_Type_free(&compress_matrix_rec_type);

    
    double  end = MPI_Wtime();

    MPI_Comm_free(&gridComm);
    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&columnComm);
    MPI_Finalize();

    if (rank == 0) {
        fprintf(stdout, "Data size: %ld x %ld x %ld\n", n1, n2, n3);
        fprintf(stdout, "Grid size: %d x %d\n", dims[Dim_X], dims[Dim_Y]);
        fprintf(stdout, "Time passed: %lf\n", end - start);
        if (correct_matrix_mult(base_A, base_B, base_C, part_mem, n1, n2, n3)) {
            fprintf(stdout, "Calculation correct\n");
        }
        else {
            fprintf(stdout, "Calculation failed\n");
        }
    }
    free(part_mem);
    free(base_mem);

    return 0;
}

#include <mpi.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <unistd.h>

double function(const double x, const double y)
{ // правая часть уравнения

    return 0;
}

double conditions(const double x, const double y)
{ //краевые условия

    if (x == 0)
        return 0;
    if (x == 1)
        return sin(y);
    if (y == 0)
        return 0;
    if (y == 1)
        return sin(x);
    else
        return 0;
    return 0;
}

double step(const size_t size)
{
    return 1.0 / (size + 1.0);
}

void first_approx_f(double **matrix, const size_t size, const double h)
{
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            matrix[i][j] = function((i + 1) * h, (j + 1) * h);
        }
    }
}

void first_approx_u(double **matrix, const size_t size, const double h)
{

    // for (size_t i = 1; i < size + 1; ++i)
    // {
    //     matrix[i][0] = conditions(i * h, 0);
    //     matrix[i][size + 1] = conditions(i * h, (size + 1) * h);
    // }
    // for (size_t j = 0; j < size + 2; ++j)
    // {
    //     matrix[0][j] = conditions(0, j * h);
    //     matrix[size + 1][j] = conditions((size + 1) * h, j * h);
    // }
    for (size_t i = 0; i < size + 2; ++i)
    {
        for (size_t j = 0; j < size + 2; ++j)
        {
            matrix[i][j] = i;
        }
    }
}

double **makeArray2D(size_t rows, size_t cols)
{
    double *data = new double[rows * cols];
    double **array2D = new double *[rows];

    for (size_t i = 0; i < rows; ++i)
        array2D[i] = &(data[i * cols]);

    return array2D;
}

void freeArray2D(double **array2D, size_t rows, size_t cols)
{
    delete[] &(array2D[0][0]);
    delete[] array2D;
}

void print2D(double **array2D, size_t rows, size_t cols)
{
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            std::cout << std::setw(2) << array2D[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[])
{
    int ProcRank, ProcSize;
    size_t N = 12;
    double h = step(N);
    double eps = 0.0001;

    int opt;
    while ((opt = getopt(argc, argv, "n:e")) != -1)
    {
        switch (opt)
        {
        case 'n':
        {
            N = std::atoi(optarg);
            break;
        }
        case 'e':
        {
            eps = std::atof(optarg);
            break;
        }
        default:
            break;
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);

    double **f = 0;
    double **u = 0;

    f = makeArray2D(N, N);
    u = makeArray2D(N + 2, N + 2);

    if (ProcRank == 0)
    {
        std::cout << "ProcSize: " << ProcSize << std::endl;

        // f = makeArray2D(N, N);
        // u = makeArray2D(N + 2, N + 2);

        first_approx_f(f, N, h);
        first_approx_u(u, N, h);
        print2D(u, N + 2, N + 2);
    }

    const int M = N / ProcSize;
    int L = (M + 2) * (N + 2);
    int *displs = new int[ProcSize];
    int *sendcounts = new int[ProcSize];

    for (size_t i = 0; i < ProcSize; i++)
    {
        displs[i] = (i == 0) ? 0 : (i * L - 2 * (N + 2));
        sendcounts[i] = L;
    }

    double **localMatrix = makeArray2D(M + 2, N + 2);

    MPI_Scatterv(&(u[0][0]), sendcounts, displs, MPI_DOUBLE, &(localMatrix[0][0]), L, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] displs;
    delete[] sendcounts;

    if (ProcRank == 1)
    {
        std::cout << "ProcRank: " << ProcRank << std::endl;
        print2D(localMatrix, M + 2, N + 2);
    }

    // Calculate



    // Clean memory
    freeArray2D(localMatrix, M + 2, N + 2);
    freeArray2D(f, N, N);
    freeArray2D(u, N + 2, N + 2);

    MPI_Finalize();

    return 0;
}

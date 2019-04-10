#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <unistd.h>

double function(const double x, const double y)
{
    return 0;
}

double conditions(const double x, const double y)
{
    // return (y - 0.5) * (y - 0.5) + (x - 0.5) * (x - 0.5);

    if (x == 0)
        return 1 - 2 * y;
    if (x == 1)
        return -1 + 2 * y;
    if (y == 0)
        return 1 - 2 * x;
    if (y == 1)
        return -1 + 2 * x;
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
    for (size_t i = 1; i < size + 1; ++i)
    {
        matrix[i][0] = conditions(i * h, 0);
        matrix[i][size + 1] = conditions(i * h, (size + 1) * h);
    }
    for (size_t j = 0; j < size + 2; ++j)
    {
        matrix[0][j] = conditions(0, j * h);
        matrix[size + 1][j] = conditions((size + 1) * h, j * h);
    }
    // for (size_t i = 0; i < size + 2; ++i)
    // {
    //     for (size_t j = 0; j < size + 2; ++j)
    //     {
    //         matrix[i][j] = i;
    //     }
    // }
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
    delete[] & (array2D[0][0]);
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
    size_t N = 8;
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

    double **locMat = makeArray2D(M + 2, N + 2);

    double u0 = 0.0, delta = 0.0;
    double dMax = 0.0;

    do
    {
        if (ProcRank == 3)
            std::cout << "0 " << ProcRank << std::endl;
        MPI_Scatterv(&(u[0][0]), sendcounts, displs, MPI_DOUBLE, &(locMat[0][0]), L, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Calculate
        if (ProcRank == 3)
            std::cout << "1 " << ProcRank << std::endl;
        for (size_t i = 1; i < M + 1; i++)
            for (size_t j = 1; j < N + 1; j++)
            {
                u0 = locMat[i][j];
                locMat[i][j] = 0.25 * (locMat[i - 1][j] + locMat[i + 1][j] + locMat[i][j - 1] + locMat[i][j + 1]); // - h * h * f[i - 1][j - 1]);
                delta = std::max(std::fabs(u0 - locMat[i][j]), delta);
            }
        if (ProcRank == 3)
            std::cout << "2 " << ProcRank << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&delta, &dMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&dMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(&(locMat[0][0]), L, MPI_DOUBLE, &(u[0][0]), sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // if (ProcRank == 0)
        // {
        //     std::cout << "D: " << dMax << std::endl;
        //     print2D(u, N + 2, N + 2);
        //     std::cout << "__________________" << std::endl;
        // }
        if (ProcRank == 3)
            std::cout << "3 " << ProcRank << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        if (ProcRank == 3)
            std::cout << "4 " << ProcRank << std::endl;
    } while (false);
    // Clean memoryp
    // delete[] displs;
    // delete[] sendcounts;
    // freeArray2D(locMat, M + 2, N + 2);
    // freeArray2D(f, N, N);
    // freeArray2D(u, N + 2, N + 2);

    MPI_Finalize();

    return 0;
}

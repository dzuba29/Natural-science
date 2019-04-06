#include <mpi.h>
#include <iostream>
#include <cmath>

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
{ //шаг для сетки

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

void first_approx_u(double **u, const size_t size, const double h)
{

    // for (size_t i = 1; i < size + 1; ++i)
    // {
    //     u[i][0] = conditions(i * h, 0);
    //     u[i][size + 1] = conditions(i * h, (size + 1) * h);
    // }
    // for (size_t j = 0; j < size + 2; ++j)
    // {
    //     u[0][j] = conditions(0, j * h);
    //     u[size + 1][j] = conditions((size + 1) * h, j * h);
    // }
    for (size_t i = 0; i < size + 2; ++i)
    {
        for (size_t j = 0; j < size + 2; ++j)
        {
            u[i][j] = 1;
        }
    }
}

// void solve(double **f, double **m, const double h, double ProcRank, double ProcSize)
// {

// }

int main(int argc, char *argv[])
{
    int ProcRank, ProcSize;
    size_t N = 10;
    double h = step(N);
    double eps = 0.0001;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);

    double *f_data = new double[N * N];
    double **f = new double *[N];
    for (size_t i = 0; i < N; ++i)
    {
        f[i] = &(f_data[i * N]);
    }

    double *u_data = new double[(N + 2) * (N + 2)];
    double **u = new double *[N + 2];
    for (size_t i = 0; i < N + 2; ++i)
    {
        u[i] = &(u_data[i * (N + 2)]);
    }

    if (ProcRank == 0)
    {
        std::cout << ProcSize << " ?????" << std::endl;
        first_approx_f(f, N, h);
        first_approx_u(u, N, h);
    }

    const int M = N / ProcSize;

    double *localMatrix_data = new double[(M + 2) * (N + 2)];
    double **localMatrix = new double *[M + 2];
    for (size_t i = 0; i < M + 2; ++i)
    {
        localMatrix[i] = &(localMatrix_data[i * (N + 2)]);
    }

    if (ProcRank == 0)
    {
        for (size_t i = 0; i < N + 2; ++i)
        {
            for (size_t j = 0; j < N + 2; ++j)
            {
                std::cout << u[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    int L = (M + 2) * (N + 2);

    // int *displs = new int[ProcSize];
    // int *sendcounts = new int[ProcSize];

    // for (size_t i = 0; i < ProcSize; i++)
    // {
    //     displs[i] = L;
    //     if (i == 0)
    //         sendcounts[i] = 0;
    //     else
    //         sendcounts[i] = i * L - 1;
    // }

    int displs[] = {0, L-2*(N + 2)};
    int sendcounts[] = {L, L};

    MPI_Scatterv(&(u[0][0]), sendcounts,displs, MPI_DOUBLE, &(localMatrix[0][0]), L, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (ProcRank == 1)
    {
        std::cout << ProcRank << " - !!!!!!!!!!!!!" << std::endl;
        for (size_t i = 0; i < M + 2; ++i)
        {
            for (size_t j = 0; j < N + 2; ++j)
            {
                std::cout << localMatrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}

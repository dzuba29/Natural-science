#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <fstream>

double function(const double x, const double y) //наш оператор Лапласа, нулевой лол
{
    return 0;
}

double conditions(const double x, const double y) //граничные условия
{
    double e = 0.0000000001;
    if (1 - e <= y && 1 + e >= y)
        return x;
    if (1 - e <= x && 1 + e >= x)
        return y * y;
    if (-e <= y && e >= y)
        return 0;
    if (-e <= x && e >= x)
        return 0;
}

double step(const size_t size) //расчет шага сетки
{
    return 1.0 / (size + 1.0);
}

void first_approx_f(double **matrix, const size_t size, const double h) //первое приближение для f
{
    for (size_t i = 0; i < size; ++i)
    {
        for (size_t j = 0; j < size; ++j)
        {
            matrix[i][j] = function((i + 1) * h, (j + 1) * h);
        }
    }
}

void first_approx_u(double **matrix, const size_t size, const double h) //первое приближение для u, заполнение граничными условиями
{
    for (size_t i = 1; i < size + 1; ++i)
    {
        matrix[i][0] = conditions(i * h, 0);
        matrix[i][size + 1] = conditions(i * h, 1);
    }
    for (size_t j = 0; j < size + 2; ++j)
    {
        matrix[0][j] = conditions(0, j * h);
        matrix[size + 1][j] = conditions(1, j * h);
    }
}

double **makeArray2D(size_t rows, size_t cols) //выделение памяти(наша матрица как бы двухмерная, но одновременно хранится как вектор)
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

void toFile(double **array2D, size_t rows, size_t cols)
{
    std::ofstream outStream("result.txt", std::ios_base::out);
    // outStream << rows << " " << cols << std::endl;
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            outStream << array2D[i][j];
            if (j != cols - 1)
                outStream << " ";
        }
        outStream << std::endl;
    }
    outStream.close();
}

int main(int argc, char *argv[])
{   
    /*
    После MPI_Init(); каждый узел выполняет то что написано в int main(), если не написано ProcRank == 0 или чему еще.
    */
    MPI_Init(&argc, &argv);
    int ProcRank, ProcSize;
    size_t N = 1000;     //размер сетки
    double h = step(N);  //шаг
    double eps = 0.0001; //точность

    int opt;
    while ((opt = getopt(argc, argv, "n:e")) != -1)
    {
        switch (opt)
        {
        case 'n':
        {
            N = std::atoi(optarg);
            h = step(N);
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

    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);

    std::ofstream DEBUG_FILE("out-" + std::to_string(ProcRank) + ".txt", std::ios_base::out);

    double **u = makeArray2D(N + 2, N + 2); //выделение памяти под апроксемирующую матрицу u

    if (ProcRank == 0)
    {
        first_approx_u(u, N, h); //заполнение начальными значениями на 0 ранге процессоров(на мастере)
    }

    const int M = N / ProcSize;                   //количество отправлямых строк матрицы u на 1 узел кластера
    int *displs_Scatterv = new int[ProcSize];     //смещение для отправляемых строк
    int *displs_Gatherv = new int[ProcSize];      //смещение для принимаемых строк
    int *sendcounts_Scatterv = new int[ProcSize]; //кол-во отправляемых строк из матрицы u каждому узлу
    int *sendcounts_Gatherv = new int[ProcSize];  // кол-во принимаемых строк с каждого узла в матриу u

    for (size_t i = 0; i < ProcSize; i++)
    {
        displs_Scatterv[i] = i * M * (N + 2);       //расчет смещения для отправки лент из матрицы u
        sendcounts_Scatterv[i] = (M + 2) * (N + 2); //расчет кол-ва отправляемых строк в ленте из матрицы u

        displs_Gatherv[i] = i * M * (N + 2) + (N + 2); //расчет смещения для принимаемых лент в матрицу u
        sendcounts_Gatherv[i] = M * (N + 2);           //расчет кол-ва принимаемых строк в ленте в матрицу u
    }
    //в sendcounts_Scatterv и sendcounts_Gatherv разные формулы т.к если оставить их одинаковыми, то при сборе результата(Gatherv) мы будем получать неверную матрицу
    double **locMat = makeArray2D(M + 2, N + 2); //выделение памяти под локальный массив в котором будут храниться отправленные строки
    double u0, lMax, dMax;                       //u0 для хранения предидущей итерации, lMax для хранения локального максимума на каждом узле, dMax - общий максимум

    do
    {
        lMax = 0.0; //локальный максимум
        /*
        Scatterv отправляет всем нодам из матрицы u со смещением displs_Scatterv и количеством sendcounts_Scatterv строки
        */

        MPI_Scatterv(&(u[0][0]), sendcounts_Scatterv, displs_Scatterv, MPI_DOUBLE, &(locMat[0][0]), sendcounts_Scatterv[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (size_t i = 1; i < M + 1; i++)
            for (size_t j = 1; j < N + 1; j++)
            {
                u0 = locMat[i][j];                                                                                 //буффер для предыдущей итерации
                locMat[i][j] = 0.25 * (locMat[i - 1][j] + locMat[i + 1][j] + locMat[i][j - 1] + locMat[i][j + 1]); //расчет по формуле
                lMax = std::max(std::fabs(u0 - locMat[i][j]), lMax);                                               //расчет локального максимума на каждом узле
            }
        /*
        Gatherv собирает все ленты с каждого узла в матрицу  u u со смещением displs_Gatherv и количеством sendcounts_Gatherv строки
        */

        MPI_Gatherv(&(locMat[0][0]) + (N + 2), sendcounts_Gatherv[0], MPI_DOUBLE, &(u[0][0]), sendcounts_Gatherv, displs_Gatherv, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Allreduce(&lMax, &dMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //Сравнивает все локальные максимумы с глобальным и если lMax > dMax, то dMax=lMax, отправляет всем узлам dMax

        MPI_Barrier(MPI_COMM_WORLD); // Синхронизация процессов чтобы они все получили и сравнили свои локальные максимумы с глобальными
    } while (dMax > eps);            //критерий останова, расчеты будут выполнятся до тех пор пока дельта между глобальным максимумом и локальным не будет превышать точность

    if (ProcRank == 0)
        toFile(u, N + 2, N + 2); // запись в файл на главном узле(мастере)

    //уборка мусора
    delete[] displs_Scatterv;
    delete[] displs_Gatherv;
    delete[] sendcounts_Scatterv;
    delete[] sendcounts_Gatherv;
    freeArray2D(locMat, M + 2, N + 2);
    freeArray2D(u, N + 2, N + 2);

    MPI_Finalize();
    return 0;
}

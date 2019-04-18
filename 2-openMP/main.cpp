#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <experimental/filesystem>
#include <unistd.h>
#include <omp.h>
#include <cmath>

double fi(double x, double y)
{
    return 0;
}

double conditions(double x, double y, double t)
{
    return pow(x, 2) - pow(y, 2);
}

double first_cond_u(double **u, size_t X, size_t Y, double t)
{
    double hX = 1.0 / X;
    double hY = 1.0 / Y;
#pragma omp parallel
    {
#pragma omp for

        for (size_t i = 1; i < Y - 1; i++)
        {
            u[i][0] = conditions(0, i * hY, t);
            u[i][X - 1] = conditions(1, i * hY, t);
        }
#pragma omp for
        for (size_t i = 0; i < X; i++)
        {
            u[0][i] = conditions(i * hX, 0, t);
            u[Y - 1][i] = conditions(i * hX, 1, t);
        }
    }
}

double first_st_u(double **u, size_t X, size_t Y)
{
    double hX = 1.0 / X;
    double hY = 1.0 / Y;
#pragma omp parallel for
    for (size_t i = 1; i < Y - 1; i++)
        for (size_t j = 1; j < X - 1; j++)
            u[i][j] = fi(j * hX, i * hY);
}

void toCSV(double **u, size_t X, size_t Y, double t, size_t N)
{
    std::ofstream outStream("res/" + std::to_string(N) + ".txt", std::ios_base::out);
    outStream << t << "\n";
    for (size_t i = 0; i < Y; ++i)
    {
        for (size_t j = 0; j < X; ++j)
        {
            outStream << u[i][j];
            if (j != Y - 1)
                outStream << " ";
        }
        outStream << std::endl;
    }
    outStream.close();
}

int main(int argc, char *argv[])
{
    size_t T = 1;   //время у.е.
    size_t nT = 10; //кол-во моментов во времени
    size_t X = 100; //размер сетки по X
    size_t Y = 100; //размер сетки по Y
    bool write_file = false;

    int opt;
    while ((opt = getopt(argc, argv, "t:n:x:y:l")) != -1)
    {
        switch (opt)
        {
        case 't':
        {
            T = std::atol(optarg);
            break;
        }
        case 'n':
        {
            nT = std::atol(optarg);
            break;
        }
        case 'x':
        {
            X = std::atol(optarg);
            break;
        }
        case 'y':
        {
            Y = std::atol(optarg);
            break;
        }
        case 'l':
        {
            write_file = true;
            break;
        }
        default:
            break;
        }
    }

    double hX = 1.0 / X; //шаг для сетки X
    double hY = 1.0 / Y; //шаг для сетки Y

    /*Выделение памяти под трехмерный массив u(t,x,y), где первая размерность указывает на время, т.е 
    в каждый момент времени t мы будем хранить расчитанную сетку X на Y, где каждая точка соотвествует температуре
    */
    double ***u = new double **[nT + 1];
    for (size_t iT = 0; iT < nT + 1; iT++)
    {
        u[iT] = new double *[Y];
        for (size_t i = 0; i < Y; i++)
            u[iT][i] = new double[X];
    }

    double t = (double)T / (double)nT; //отрезок времени
    double l1 = t / (hX * hX);         //коэф лямбда1 из формулы
    double l2 = t / (hY * hY);         //коэф лямбда2 из формулы
    double c = 1 - 2 * (l1 + l2);      //сократили запись(формула)

    // Условие устойчивости
    if (l1 + l2 > 0.5)
    {
        std::cout << "Stability condition not satisfied!\n";
        return 0;
    }

    auto begTime = std::chrono::steady_clock::now();


    first_st_u(u[0], X, Y);      //заполнение в 0 момент времени начальных условий
    first_cond_u(u[0], X, Y, 0); //заполнение граничных условий в 0 момент времени

    for (size_t iT = 1; iT < nT + 1; iT++)
    {
#pragma omp parallel for
        for (size_t i = 1; i < Y - 1; i++)
            for (size_t j = 1; j < X - 1; j++)
                /*
                Расчет по формуле 
                https://3ys.ru/metody-resheniya-differentsialnykh-uravnenij/uravnenie-teploprovodnosti.html
                под 2.78
                в iT момент времени 
                */
                u[iT][i][j] = c * u[iT - 1][i][j] + l1 * (u[iT - 1][i][j + 1] + u[iT - 1][i][j - 1]) + l2 * (u[iT - 1][i + 1][j] + u[iT - 1][i - 1][j]);

        // Пересчет граничных условий в iT момент времени
        first_cond_u(u[iT], X, Y, t * iT);
    }

    double calc_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - begTime).count();

    std::cout << omp_get_max_threads() << " " << calc_time << std::endl;

    // вывод решения в файл
    if (write_file)
    {
        std::experimental::filesystem::remove_all("res");
        std::experimental::filesystem::create_directory("res");
        for (size_t iT = 0; iT < nT + 1; iT++)
        {
            toCSV(u[iT], X, Y, t * iT, iT);
        }
    }

    // Освобождение памяти
    for (size_t iT = 0; iT < nT + 1; iT++)
    {
        for (size_t i = 0; i < Y; i++)
            delete[] u[iT][i];
        delete[] u[iT];
    }
    delete[] u;

    return 0;
}

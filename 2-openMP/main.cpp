#include <iostream>
#include <fstream>
#include <string>

double fi(double x, double y)
{
    return x * y * 0.2;
}

double conditions(double x, double y, double t)
{
    double e = 0.000001;
    if (1 - e <= y && 1 + e >= y)
        return x;
    if (1 - e <= x && 1 + e >= x)
        return y * y;
    if (-e <= y && e >= y)
        return 0;
    if (-e <= x && e >= x)
        return 0;
}

double fillCondi(double **u, size_t X, size_t Y, double t)
{
    double hX = 1.0 / (X + 1.0);
    double hY = 1.0 / (Y + 1.0);
    for (size_t i = 1; i < Y - 1; i++)
    {
        u[i][0] = conditions(0, i * hY, t);
        u[i][X - 1] = conditions(1, i * hY, t);
    }
    for (size_t i = 0; i < X; i++)
    {
        u[0][i] = conditions(i * hX, 0, t);
        u[Y - 1][i] = conditions(i * hX, 1, t);
    }
}

void makeCSV(double **u, size_t X, size_t Y, double t, size_t N)
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

int main()
{
    double T = 1;
    size_t nT = 10;

    size_t X = 500;
    size_t Y = 500;

    double hX = 1.0 / (X + 1.0);
    double hY = 1.0 / (Y + 1.0);

    // Выделение памяти
    double ***cube = new double **[nT + 1];
    for (size_t iT = 0; iT < nT + 1; iT++)
    {
        cube[iT] = new double *[Y];
        for (size_t i = 0; i < Y; i++)
            cube[iT][i] = new double[X];
    }

    // Начальное условие
    for (size_t i = 1; i < Y - 1; i++)
        for (size_t j = 1; j < X - 1; j++)
            cube[0][i][j] = fi(j * hX, i * hY);
    fillCondi(cube[0], X, Y, 0);

    double l1 = T / nT / hX * hX;
    double l2 = T / nT / hY * hY;
    double c = 1 - 2 * (l1 + l2);
    double r = T / nT;

    makeCSV(cube[0], X, Y, 0, 0);

    for (size_t iT = 1; iT < nT + 1; iT++)
    {   
        #pragma omp parallel for
        for (size_t i = 1; i < Y - 1; i++)
            for (size_t j = 1; j < X - 1; j++)
                cube[iT][i][j] = c * cube[iT - 1][i][j] + l1 * (cube[iT - 1][i][j + 1] + cube[iT - 1][i][j - 1]) + l2 * (cube[iT - 1][i + 1][j] + cube[iT - 1][i - 1][j]);

        // заполнение границ
        fillCondi(cube[iT], X, Y, r * iT);
        makeCSV(cube[iT], X, Y, r * iT, iT);
    }

    // Освлобождение памяти
    for (size_t iT = 0; iT < nT + 1; iT++)
    {
        for (size_t i = 0; i < Y; i++)
            delete[] cube[iT][i];
        delete[] cube[iT];
    }
    delete[] cube;

    return 0;
}

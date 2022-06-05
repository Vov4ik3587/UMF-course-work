#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

//ТУТ ПОЯСНИТЕЛЬНАЯ ЗАПИСКА//
/*В файле input.txt задаётся расчётная область: граничные точки по х, граничные точки по у, количество узлов по х, количество узлов по у, коэффициент разрядки по х, коэффициент разрядки по y
В файле timeGrid.txt отрезок времени: граничные точки, количество узлов, коэффициент разрядки.
  В файлах firstboundary.txt, secondboundary.txt, thirdboundary.txt находятся номера границ, на которых заданы соответствующие краевые условия
   Для учёта краевых условий пронумеровали от 0 до 3 границы расчётной области против часовой стрелки начиная с левой границы
   Функции и коэффициенты для краевых условий задаются внутри программы в соответствующих функциях для нужных границ
   Формат хранения глобальной матрицы – разреженный
   Для решения используется Локально-оптимальная схема с предобуславливание LU(sq)*/

class slae
{
private:
    int _nx;              // количество узлов по оси х
    int _ny;              // количество узлов по оси у
    std::vector<int> _ig; // Профиль глобальной матрицы
    std::vector<int> _jg;
    std::vector<double> _al;                                      //Нижний треугольник матрицы
    std::vector<double> _au;                                      //Верхний треугольник матрицы
    std::vector<double> _di;                                      //Диагональные элементы матрицы
    std::vector<std::vector<double>> _loc_mass;                   //Локальная матрица масс
    std::vector<std::vector<std::vector<double>>> _loc_stiffness; // Локальная матрицы жёсткости
    std::vector<std::vector<double>> _grid;                       //Пронумерованные точки в узлах сетки
    std::vector<double> _time;                                    //Сетка по времени
    std::vector<std::vector<std::vector<int>>> _el_bonds;         // Связи узлов и элементов
    std::vector<double> _glob_b;                                  //Вектор правой части
    std::vector<std::vector<double>> _q;                          //Вектор решения СЛАУ (коэффициенты разложения искомой функции по биквадратичному базису)
public:
    slae()
    {
        this->_nx = 0;
        this->_ny = 0;
        this->_ig = std::vector<int>();
        this->_jg = std::vector<int>();
        this->_al = std::vector<double>();
        this->_au = std::vector<double>();
        this->_di = std::vector<double>();
        this->_loc_mass = std::vector<std::vector<double>>();
        this->_loc_stiffness = std::vector<std::vector<std::vector<double>>>();
        this->_grid = std::vector<std::vector<double>>();
        this->_el_bonds = std::vector<std::vector<std::vector<int>>>();
        this->_glob_b = std::vector<double>();
        this->_q = std::vector<std::vector<double>>();
        this->_time = std::vector<double>();
    }

    //Вовзращает значение функции правой части уравнение в точке
    double f(double x, double y, double t)
    {
        return 0;
    }

    //Возвращает значение функции коэффициента диффузии в точке
    double lambda(double x, double y)
    {
        return 1.;
    }

    //Возвращает значение коэффициента гамма/сигма/хи
    double gamma(double x, double y)
    {
        return 0.;
    }

    double sigma(double x, double y)
    {
        return 1.;
    }

    double khi(double x, double y)
    {
        return 1.;
    }

    //Возвращает значение функции Тета во втором краевом условии
    double theta(double x, double y, int Border, double t)
    {
        switch (Border)
        {
        case 0:
        {
            return -3.0;
            break;
        }
        case 1:
        {
            return -2.;
            break;
        }
        case 2:
        {
            return 4.0;
            break;
        }
        case 3:
        {
            return 0.0;
            break;
        }
        }
    }

    //Возвращает значение коэффициента бета в третьем краевом условии
    double beta()
    {
        return 0.5;
    }

    //Возвращает значение функции Убета В третьем краевом условии
    double Ubeta(double x, double y, int Border, double t)
    {
        switch (Border)
        {
        case 0:
        {
            return x * x - 4.;
            break;
        }
        case 1:
        {
            return 0.0;
            break;
        }
        case 2:
        {
            return 3.0;
            break;
        }
        case 3:
        {
            return x;
            break;
        }
        }
    }

    //Возвращает значение функции, задающей первое краевое условие
    double Ug(double x, double y, int Border, double t)
    {
        switch (Border)
        {
        case 0:
        {
            return t + y * y * y * y * y;
            break;
        }
        case 1:
        {
            return t + x * x * x * x * x;
            break;
        }
        case 2:
        {
            return t + 1. + y * y * y * y * y;
            break;
        }
        case 3:
        {
            return t + 1. + x * x * x * x * x;
            break;
        }
        }
    }

    //Возвращает значение искомой функции в точке для проверки точности решения
    double U(double x, double y, double t)
    {

        return x * x * x * x * x + y * y * y * y * y + t;
    }

    //Начальное значение
    double U_0(double x, double y)
    {
        return x * x;
    }

    //Начальная производная
    double dU_0(double x, double y)
    {
        return 0.;
    }
    //Значение функции на первом шаге
    double U_1(double x, double y)
    {
        return x * x;
    }

    void set_profile(std::vector<std::vector<int>> NodeBonds)
    {
        _ig.resize((2 * _nx - 1) * (2 * _ny - 1) + 1);
        _di.resize((2 * _nx - 1) * (2 * _ny - 1));
        _ig[0] = 0;
        _ig[1] = 0;
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            int k = 0;
            for (int j = 0; j < NodeBonds[i].size(); j++)
            {
                k++;
                _jg.push_back(NodeBonds[i][j]);
            }
            _ig[i + 1] = _ig[i] + k;
        }
    }

    void making_time_grid()
    {
        std::fstream input;
        input.open("timeGrid.txt");
        double tmin, tmax;
        int nt;
        input >> tmin >> tmax >> nt;
        double ht = (tmax - tmin) / (nt - 1);
        _time.resize(nt);
        for (int i = 0; i < nt; i++)
        {
            _time[i] = tmin + ht * i;
        }
        input.close();
        _q.resize(nt);
    }
    //Создаёт сетку и список связи элементов и узлов
    void making_grid()
    {
        std::fstream input;
        input.open("input.txt");
        double xmin, xmax, ymin, ymax;
        double kx, ky;
        input >> xmin;
        input >> xmax;
        input >> ymin;
        input >> ymax;
        input >> _nx;
        input >> _ny;
        input >> kx;
        input >> ky;
        std::vector<double> x;
        std::vector<double> y;
        if (kx == 1)
        {
            double hx = (xmax - xmin) / (_nx - 1);
            x.resize(_nx);
            for (int i = 0; i < _nx; i++)
            {
                x[i] = xmin + i * hx;
            }
        }
        else
        {
            double hx = (xmax - xmin) / ((1. - pow(kx, _nx - 1)) / (1. - kx));
            x.resize(_nx);
            for (int i = 0; i < _nx; i++)
            {
                x[i] = xmin + hx * ((1. - pow(kx, i)) / (1. - kx));
            }
        }

        if (ky == 1)
        {
            double hy = (ymax - ymin) / (_ny - 1);
            y.resize(_ny);
            for (int i = 0; i < _ny; i++)
            {
                y[i] = ymin + i * hy;
            }
        }
        else
        {
            double hy = (ymax - ymin) / ((1. - pow(ky, _ny - 1)) / (1. - ky));
            y.resize(_ny);
            for (int i = 0; i < _ny; i++)
            {
                y[i] = ymin + hy * ((1. - pow(ky, i)) / (1. - ky));
            }
        }

        _grid.resize(_nx * _ny);
        int k = 0;
        for (int i = 0; i < _ny; i++)
            for (int j = 0; j < _nx; j++)
            {
                _grid[k].resize(2);
                _grid[k][0] = x[j];
                _grid[k][1] = y[i];
                k++;
            }
        input.close();
        _el_bonds.resize(2);
        _el_bonds[0].resize((_nx - 1) * (_ny - 1));
        _el_bonds[1].resize((_nx - 1) * (_ny - 1));
        for (int i = 0; i < (_ny - 1); i++)
        {
            for (int j = 0; j < (_nx - 1); j++)
            {
                _el_bonds[0][i * (_nx - 1) + j].resize(9);
                _el_bonds[0][i * (_nx - 1) + j][0] = 2 * i * (2 * _nx - 1) + j * 2;
                _el_bonds[0][i * (_nx - 1) + j][1] = 2 * i * (2 * _nx - 1) + j * 2 + 1;
                _el_bonds[0][i * (_nx - 1) + j][2] = 2 * i * (2 * _nx - 1) + j * 2 + 2;
                _el_bonds[0][i * (_nx - 1) + j][3] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx - 1;
                _el_bonds[0][i * (_nx - 1) + j][4] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx;
                _el_bonds[0][i * (_nx - 1) + j][5] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx + 1;
                _el_bonds[0][i * (_nx - 1) + j][6] = 2 * (i + 1) * (2 * _nx - 1) + j * 2;
                _el_bonds[0][i * (_nx - 1) + j][7] = 2 * (i + 1) * (2 * _nx - 1) + j * 2 + 1;
                _el_bonds[0][i * (_nx - 1) + j][8] = 2 * (i + 1) * (2 * _nx - 1) + j * 2 + 2;
                _el_bonds[1][i * (_nx - 1) + j].resize(4);
                _el_bonds[1][i * (_nx - 1) + j][0] = i * _nx + j;
                _el_bonds[1][i * (_nx - 1) + j][1] = i * _nx + j + 1;
                _el_bonds[1][i * (_nx - 1) + j][2] = (i + 1) * _nx + j;
                _el_bonds[1][i * (_nx - 1) + j][3] = (i + 1) * _nx + j + 1;
            }
        }

        std::vector<std::vector<int>> NodeBonds;
        NodeBonds.resize((2 * _nx - 1) * (2 * _ny - 1));
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            for (int j = 0; j < 9; j++)
            {
                for (int k = 0; k < 9; k++)
                {
                    if (_el_bonds[0][i][j] > _el_bonds[0][i][k])
                    {
                        NodeBonds[_el_bonds[0][i][j]].push_back(_el_bonds[0][i][k]);
                    }
                }
            }
        }

        for (int i = 0; i < ((2 * _nx - 1) * (2 * _ny - 1)); i++)
        {
            sort(NodeBonds[i].begin(), NodeBonds[i].end());
            auto last = std::unique(NodeBonds[i].begin(), NodeBonds[i].end());
            NodeBonds[i].erase(last, NodeBonds[i].end());
        }
        set_profile(NodeBonds);
    }

    void init_hyperbolic_approx(int j)
    {
        _q[j].resize((2 * _nx - 1) * (2 * _ny - 1));
        std::fill(_q[j].begin(), _q[j].end(), 1.);
    }

    void making_initial_conds()
    {
        double xmin, ymin, hx, hy;
        _q[0].resize((2 * _nx - 1) * (2 * _ny - 1));
        _q[1].resize((2 * _nx - 1) * (2 * _ny - 1));
        xmin = _grid[0][0];
        ymin = _grid[0][1];
        hx = (_grid[1][0] - _grid[0][0]) / 2.;
        hy = (_grid[_nx][1] - _grid[0][1]) / 2.;
        int k = 0;
        for (int i = 0; i < 2 * _ny - 1; i++)
        {
            for (int j = 0; j < 2 * _nx - 1; j++)
            {
                _q[0][k] = U_0(xmin + hx * j, ymin + hy * i);
                _q[1][k] = _q[0][k] + dU_0(xmin + hx * j, ymin + hy * i) * (_time[1] - _time[0]);
                //_Q[1][k] = U_1(xmin + hx * j, ymin + hy * i);
                k++;
            }
        }
    }

    //Производит сборку локальных матриц
    void local_build(double hx, double hy)
    {
        std::vector<std::vector<std::vector<double>>> LocMassPsi = std::vector<std::vector<std::vector<double>>>(3);
        std::vector<std::vector<std::vector<double>>> LocGPsi = std::vector<std::vector<std::vector<double>>>(3);
        for (int i = 0; i < 3; i++)
        {
            LocMassPsi[i].resize(3);
        }
        //Матрица для x/y
        LocMassPsi[0][0].resize(3);
        LocMassPsi[0][1].resize(3);
        LocMassPsi[0][2].resize(3);
        LocMassPsi[0][0][0] = 1.0 / 60.0;
        LocMassPsi[0][0][1] = 0.;
        LocMassPsi[0][0][2] = -1.0 / 60.0;
        LocMassPsi[0][1][0] = 0.;
        LocMassPsi[0][1][1] = 4.0 / 15.0;
        LocMassPsi[0][1][2] = 1.0 / 15.0;
        LocMassPsi[0][2][0] = -1.0 / 60.0;
        LocMassPsi[0][2][1] = 1.0 / 15.0;
        LocMassPsi[0][2][2] = 7.0 / 60.0;

        //Матрица для 1-x/1-y
        LocMassPsi[1][0].resize(3);
        LocMassPsi[1][1].resize(3);
        LocMassPsi[1][2].resize(3);
        LocMassPsi[1][0][0] = 7.0 / 60.0;
        LocMassPsi[1][0][1] = 1.0 / 15.0;
        LocMassPsi[1][0][2] = -1.0 / 60.0;
        LocMassPsi[1][1][0] = 1.0 / 15.0;
        LocMassPsi[1][1][1] = 4.0 / 15.0;
        LocMassPsi[1][1][2] = 0.;
        LocMassPsi[1][2][0] = -1.0 / 60.0;
        LocMassPsi[1][2][1] = 0.;
        LocMassPsi[1][2][2] = 1.0 / 60.0;

        //Стандартная матрица
        LocMassPsi[2][0].resize(3);
        LocMassPsi[2][1].resize(3);
        LocMassPsi[2][2].resize(3);
        LocMassPsi[2][0][0] = 2.0 / 15.0;
        LocMassPsi[2][0][1] = 1.0 / 15.0;
        LocMassPsi[2][0][2] = -1.0 / 30.0;
        LocMassPsi[2][1][0] = 1.0 / 15.0;
        LocMassPsi[2][1][1] = 8.0 / 15.0;
        LocMassPsi[2][1][2] = 1.0 / 15.0;
        LocMassPsi[2][2][0] = -1.0 / 30.0;
        LocMassPsi[2][2][1] = 1.0 / 15.0;
        LocMassPsi[2][2][2] = 2.0 / 15.0;

        for (int i = 0; i < 3; i++)
        {
            LocGPsi[i].resize(3);
        }

        //Матрица для x/y
        LocGPsi[0][0].resize(3);
        LocGPsi[0][1].resize(3);
        LocGPsi[0][2].resize(3);
        LocGPsi[0][0][0] = 1. / 2.;
        LocGPsi[0][0][1] = -2.0 / 3.;
        LocGPsi[0][0][2] = 1.0 / 6.0;
        LocGPsi[0][1][0] = -2.0 / 3.0;
        LocGPsi[0][1][1] = 8.0 / 3.0;
        LocGPsi[0][1][2] = -2.;
        LocGPsi[0][2][0] = 1. / 6.;
        LocGPsi[0][2][1] = -2.;
        LocGPsi[0][2][2] = 11. / 6.;

        //Матрица для 1-x/1-y
        LocGPsi[1][0].resize(3);
        LocGPsi[1][1].resize(3);
        LocGPsi[1][2].resize(3);
        LocGPsi[1][0][0] = 11.0 / 6.0;
        LocGPsi[1][0][1] = -2.;
        LocGPsi[1][0][2] = 1.0 / 6.0;
        LocGPsi[1][1][0] = -2.;
        LocGPsi[1][1][1] = 8.0 / 3.0;
        LocGPsi[1][1][2] = -2. / 3.;
        LocGPsi[1][2][0] = 1. / 6.;
        LocGPsi[1][2][1] = -2. / 3.;
        LocGPsi[1][2][2] = 1. / 2.;

        LocGPsi[2][0].resize(3);
        LocGPsi[2][1].resize(3);
        LocGPsi[2][2].resize(3);
        LocGPsi[2][0][0] = 7. / 3.;
        LocGPsi[2][0][1] = -8. / 3.;
        LocGPsi[2][0][2] = 1. / 3.;
        LocGPsi[2][1][0] = -8. / 3.;
        LocGPsi[2][1][1] = 16. / 3.;
        LocGPsi[2][1][2] = -8. / 3.;
        LocGPsi[2][2][0] = 1. / 3.;
        LocGPsi[2][2][1] = -8. / 3.;
        LocGPsi[2][2][2] = 7. / 3.;

        _loc_mass.resize(9);
        for (int i = 0; i < 9; i++)
        {
            _loc_mass[i].resize(9);
            for (int j = 0; j < 9; j++)
            {
                _loc_mass[i][j] = hx * hy * LocMassPsi[2][i / 3][j / 3] * LocMassPsi[2][i % 3][j % 3];
            }
        }
        _loc_stiffness.resize(5);
        for (int k = 0; k < 5; k++)
        {
            _loc_stiffness[k].resize(9);
        }
        for (int i = 0; i < 9; i++)
        {
            _loc_stiffness[0][i].resize(9);
            _loc_stiffness[1][i].resize(9);
            _loc_stiffness[2][i].resize(9);
            _loc_stiffness[3][i].resize(9);
            _loc_stiffness[4][i].resize(9);
            for (int j = 0; j < 9; j++)
            {
                _loc_stiffness[0][i][j] = hy * LocGPsi[0][i % 3][j % 3] * LocMassPsi[0][i / 3][j / 3] / hx + hx * LocGPsi[0][i / 3][j / 3] * LocMassPsi[0][i % 3][j % 3] / hy;
                _loc_stiffness[1][i][j] = hy * LocGPsi[1][i % 3][j % 3] * LocMassPsi[0][i / 3][j / 3] / hx + hx * LocGPsi[0][i / 3][j / 3] * LocMassPsi[1][i % 3][j % 3] / hy;
                _loc_stiffness[2][i][j] = hy * LocGPsi[0][i % 3][j % 3] * LocMassPsi[1][i / 3][j / 3] / hx + hx * LocGPsi[1][i / 3][j / 3] * LocMassPsi[0][i % 3][j % 3] / hy;
                _loc_stiffness[3][i][j] = hy * LocGPsi[1][i % 3][j % 3] * LocMassPsi[1][i / 3][j / 3] / hx + hx * LocGPsi[1][i / 3][j / 3] * LocMassPsi[1][i % 3][j % 3] / hy;
                _loc_stiffness[4][i][j] = hy * LocGPsi[2][i % 3][j % 3] * LocMassPsi[2][i / 3][j / 3] / hx + hx * LocGPsi[2][i / 3][j / 3] * LocMassPsi[2][i % 3][j % 3] / hy;
            }
        }
    }

    //Производит сборку глобальной матрицы по элементам
    void elliptic_global_build(double t)
    {
        double gamma = 0.;
        double hx, hy;
        std::vector<double> LocB = std::vector<double>(9);
        std::vector<std::vector<double>> LocG;
        _glob_b.resize((2 * _nx - 1) * (2 * _ny - 1));
        _al.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        _au.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            hx = _grid[_el_bonds[1][i][1]][0] - _grid[_el_bonds[1][i][0]][0];
            hy = _grid[_el_bonds[1][i][2]][1] - _grid[_el_bonds[1][i][0]][1];
            local_build(hx, hy);
            get_local_gamma(i, gamma);
            get_local_f(LocB, i, t);
            mass_mx_mult_vec(LocB);
            for (int k = 0; k < 9; k++)
            {
                _glob_b[_el_bonds[0][i][k]] += LocB[k];
            }
            get_loc_G(i, LocG);
            for (int k = 0; k < 9; k++)
            {
                _di[_el_bonds[0][i][k]] += gamma * _loc_mass[k][k] + LocG[k][k];
            }

            int Index;
            for (int k = 1; k < 9; k++)
            {
                for (int j = 0; j < k; j++)
                {
                    get_index(_el_bonds[0][i][k], _el_bonds[0][i][j], Index);
                    _al[Index] += gamma * _loc_mass[k][j] + LocG[k][j];
                    _au[Index] += gamma * _loc_mass[j][k] + LocG[j][k];
                }
            }
        }
    }

    //Глобальная сборка для трёхслойной схемы
    void three_layer_global_bild(double t0, double t1, double t2)
    {
        double gamma = 0.;
        double sigma = 0.;
        double khi = 0.;
        double hx, hy;

        std::vector<double> loc_b = std::vector<double>(9);
        std::vector<double> q_0 = std::vector<double>(9);
        std::vector<double> q_1 = std::vector<double>(9);
        std::vector<double> q_1_G = std::vector<double>(9);
        std::vector<std::vector<double>> LocG;

        double dt = t2 - t0;
        double dt0 = t2 - t1;
        double dt1 = t1 - t0;

        _glob_b.resize((2 * _nx - 1) * (2 * _ny - 1));
        _al.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        _au.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            hx = _grid[_el_bonds[1][i][1]][0] - _grid[_el_bonds[1][i][0]][0];
            hy = _grid[_el_bonds[1][i][2]][1] - _grid[_el_bonds[1][i][0]][1];

            local_build(hx, hy);

            get_local_gamma(i, gamma);
            get_local_sigma(i, sigma);
            get_local_khi(i, khi);
            get_local_f(loc_b, i, t1);

            for (int k = 0; k < 9; k++)
            {
                q_0[k] = _q[0][_el_bonds[0][i][k]];
                q_1[k] = _q[1][_el_bonds[0][i][k]];
            }

            q_1_G = q_1;
            mass_mx_mult_vec(loc_b);
            mass_mx_mult_vec(q_0);
            mass_mx_mult_vec(q_1);
            get_loc_G(i, LocG);
            stiffness_mx_mult_vec(q_1_G, LocG);

            for (int k = 0; k < 9; k++)
            {
                _glob_b[_el_bonds[0][i][k]] += loc_b[k] - khi * (2. * q_0[k] / (dt * dt1) - 2. * q_1[k] / (dt1 * dt0)) -
                                               sigma * ((dt0 - dt1) * q_1[k] / (dt1 * dt0) - dt0 * q_0[k] / (dt * dt1)) -
                                               gamma * q_1[k] - q_1_G[k];
            }

            for (int k = 0; k < 9; k++)
            {
                _di[_el_bonds[0][i][k]] += 2. * khi * _loc_mass[k][k] / (dt * dt0) + sigma * dt1 * _loc_mass[k][k] / (dt0 * dt);
            }

            int Index;
            for (int k = 1; k < 9; k++)
            {
                for (int j = 0; j < k; j++)
                {
                    get_index(_el_bonds[0][i][k], _el_bonds[0][i][j], Index);
                    _al[Index] += 2. * khi * _loc_mass[k][j] / (dt * dt0) + sigma * dt1 * _loc_mass[k][j] / (dt0 * dt);
                    _au[Index] += 2. * khi * _loc_mass[k][j] / (dt * dt0) + sigma * dt1 * _loc_mass[k][j] / (dt0 * dt);
                }
            }
        }
    }

    void get_glob_b(int j)
    {
        double hx, hy, x, y;
        int size = _el_bonds[0].size();
        std::vector<std::vector<double>> fictive_grid = std::vector<std::vector<double>>((2 * _nx - 1) * (2 * _ny - 1));
        for (int i = 0; i < size; i++)
        {
            x = _grid[_el_bonds[1][i][0]][0];
            y = _grid[_el_bonds[1][i][0]][1];
            hx = (_grid[_el_bonds[1][i][1]][0] - _grid[_el_bonds[1][i][0]][0]) / 2.;
            hy = (_grid[_el_bonds[1][i][2]][1] - _grid[_el_bonds[1][i][0]][1]) / 2.;
            for (int k = 0; k < 9; k++)
            {
                fictive_grid[_el_bonds[0][i][k]].resize(2);
                fictive_grid[_el_bonds[0][i][k]][0] = x + hx * (k % 3);
                fictive_grid[_el_bonds[0][i][k]][1] = y + hy * (k / 3);
            }
        }
        size = fictive_grid.size();
        for (int i = 0; i < size; i++)
        {
            _glob_b[i] = f(fictive_grid[i][0], fictive_grid[i][1], _time[j - 1]);
        }
    }

    //умножение матрицы на вектор
    void mx_mult_vec(std::vector<double> al, std::vector<double> au, std::vector<double> di, std::vector<double> &f)
    {
        std::vector<double> res(f.size());
        for (int i = 0; i < f.size(); i++)
        {
            res[i] = di[i] * f[i];
            for (int j = _ig[i]; j < _ig[i + 1]; j++)
            {
                res[i] += al[j] * f[_jg[j]];
                res[_jg[j]] += au[j] * f[i];
            }
        }
        f = res;
    }

    //Вывод матрицы для проверки
    void print_mass_mx()
    {
        for (int i = 0; i < 9; i++)
        {
            for (int j = 0; j < 9; j++)
            {
                std::cout << _loc_mass[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    //Умножение локальной матрицы жёсткости на вектор
    void stiffness_mx_mult_vec(std::vector<double> &vec, std::vector<std::vector<double>> LocG)
    {
        std::vector<double> res = std::vector<double>(9);
        for (int i = 0; i < 9; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < 9; j++)
            {
                sum += vec[i] * LocG[i][j];
            }
            res[i] = sum;
        }
        vec = res;
    }

    //Вовзращает положение элемента в векторах al/au для его вставки в глобальную матрицу
    void get_index(int i, int j, int &Index)
    {
        Index = _ig[i];
        while (_jg[Index] != j)
        {
            Index++;
        }
    }

    //Собирает из локальных матриц жёсткости, полученных в ходе разложение коэффициента диффузии по билинейным базисным функциям
    void get_loc_G(int num, std::vector<std::vector<double>> &LocG)
    {
        std::vector<double> LocPhiLam = std::vector<double>(4);
        get_loc_lam(num, LocPhiLam);
        LocG.resize(9);
        for (int i = 0; i < 9; i++)
        {
            LocG[i].resize(9);
            for (int j = 0; j < 9; j++)
            {
                LocG[i][j] = LocPhiLam[0] * _loc_stiffness[0][i][j] + LocPhiLam[1] * _loc_stiffness[1][i][j] + LocPhiLam[2] * _loc_stiffness[2][i][j] + LocPhiLam[3] * _loc_stiffness[3][i][j];
            }
        }

        /*double lam = 0.;
        for (int i = 0; i < 4; i++)
        {
           lam += LocPhiLam[i];
        }
        lam /= 4;

        for (int i = 0; i < 9; i++)
        {
           LocG[i].resize(9);
           for (int j = 0; j < 9; j++)
           {
              LocG[i][j] = lam * _LocG[4][i][j];
           }
        }*/
    }

    //Возвращает коэффициенты разложения коэффициента диффузии по билинейным базисным функциям
    void get_loc_lam(int num, std::vector<double> &LocPhiLam)
    {
        for (int i = 0; i < 4; i++)
        {
            LocPhiLam[i] = lambda(_grid[_el_bonds[1][num][i]][0], _grid[_el_bonds[1][num][i]][1]);
        }
        int k = 0;
    }

    //Возвращают усреднённые коэффициенты при матрице масс для элемента
    void get_local_gamma(int num, double &gam)
    {
        double hx, hy;
        double g = 0.0;
        hx = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.;
        hy = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            g += gamma(_grid[_el_bonds[1][num][0]][0] + hx * (i % 3), _grid[_el_bonds[1][num][0]][1] + hy * (i / 3));
        }
        gam = g / 9.;
    }

    void get_local_sigma(int num, double &sig)
    {
        double hx, hy;
        double s = 0.0;
        hx = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.;
        hy = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            s += sigma(_grid[_el_bonds[1][num][0]][0] + hx * (i % 3), _grid[_el_bonds[1][num][0]][1] + hy * (i / 3));
        }
        sig = s / 9.;
    }

    void get_local_khi(int num, double &kh)
    {
        double hx, hy;
        double k = 0.0;
        hx = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.;
        hy = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            k += khi(_grid[_el_bonds[1][num][0]][0] + hx * (i % 3), _grid[_el_bonds[1][num][0]][1] + hy * (i / 3));
        }
        kh = k / 9.;
    }

    //Умножает вектор на матрицу масс для вставки в правую часть СЛАУ
    void mass_mx_mult_vec(std::vector<double> &f)
    {
        std::vector<double> res = std::vector<double>(9);
        for (int i = 0; i < 9; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < 9; j++)
            {
                sum += f[i] * _loc_mass[i][j];
            }
            res[i] = sum;
        }
        f = res;
    }

    //Получает вектор локальных значений правой части на конечном элементе
    void get_local_f(std::vector<double> &LocB, int num, double t)
    {
        double hx, hy;
        hx = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.;
        hy = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.;
        std::vector<double> res = std::vector<double>(9);
        for (int i = 0; i < 9; i++)
        {
            double x = _grid[_el_bonds[1][num][0]][0] + hx * (i % 3);
            double y = _grid[_el_bonds[1][num][0]][1] + hy * (i / 3);
            res[i] = f(x, y, t);
        }
        LocB = res;
    }

    //Учёт второго краевого условия
    void second_boundary(double t)
    {
        std::fstream input;
        input.open("secondboundary.txt");
        std::vector<int> BoundaryBorder;
        int variable;
        while (input >> variable)
        {
            BoundaryBorder.push_back(variable);
        }
        input.close();
        for (int Border : BoundaryBorder)
        {
            switch (Border)
            {
            case 0:
            {
                std::vector<double> LocTheta = std::vector<double>(3);
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1);
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]);
                    get_local_theta(num, h / 2.0, Border, LocTheta, t);
                    _glob_b[_el_bonds[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][3]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][6]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
                }
                break;
            }
            case 1:
            {
                std::vector<double> LocTheta = std::vector<double>(3);
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i;
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _glob_b[_el_bonds[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][1]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][2]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
                }
                break;
            }
            case 2:
            {
                std::vector<double> LocTheta = std::vector<double>(3);
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1) + _ny - 2;
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _glob_b[_el_bonds[0][num][2]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][5]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
                }
                break;
            }
            case 3:
            {
                std::vector<double> LocTheta = std::vector<double>(3);
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i + (_nx - 1) * (_ny - 2);
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _glob_b[_el_bonds[0][num][6]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][7]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
                }
                break;
            }
            }
        }
    }

    //Возвращает локальный вектор значений функции тета на ребре при учёте второго краевого условия
    void get_local_theta(int num, double h, int Border, std::vector<double> &LocTheta, double t)
    {

        for (int i = 0; i < 3; i++)
        {
            switch (Border)
            {
            case 0:
            {
                LocTheta[i] = theta(_grid[_el_bonds[1][num][0]][0], _grid[_el_bonds[1][num][0]][1] + h * i, Border, t);
                break;
            }
            case 1:
            {
                LocTheta[i] = theta(_grid[_el_bonds[1][num][0]][0] + h * i, _grid[_el_bonds[1][num][0]][1], Border, t);
                break;
            }
            case 2:
            {
                LocTheta[i] = theta(_grid[_el_bonds[1][num][1]][0], _grid[_el_bonds[1][num][1]][1] + h * i, Border, t);
                break;
            }
            case 3:
            {
                LocTheta[i] = theta(_grid[_el_bonds[1][num][2]][0] + h * i, _grid[_el_bonds[1][num][2]][1], Border, t);
                break;
            }
            default:
                break;
            }
        }
    }

    //Возвращает локальный вектор значений функции Убета на ребре для учёта в третем краевом условии
    void get_local_Ubeta(int num, double h, int Border, std::vector<double> &LocUbeta, double t)
    {
        for (int i = 0; i < 3; i++)
        {
            switch (Border)
            {
            case 0:
            {
                LocUbeta[i] = Ubeta(_grid[_el_bonds[1][num][0]][0], _grid[_el_bonds[1][num][0]][1] + h * i, Border, t);
                break;
            }
            case 1:
            {
                LocUbeta[i] = Ubeta(_grid[_el_bonds[1][num][0]][0] + h * i, _grid[_el_bonds[1][num][0]][1], Border, t);
                break;
            }
            case 2:
            {
                LocUbeta[i] = Ubeta(_grid[_el_bonds[1][num][1]][0], _grid[_el_bonds[1][num][1]][1] + h * i, Border, t);
                break;
            }
            case 3:
            {
                LocUbeta[i] = Ubeta(_grid[_el_bonds[1][num][2]][0] + h * i, _grid[_el_bonds[1][num][2]][1], Border, t);
                break;
            }
            default:
                break;
            }
        }
    }

    //Учёт третьего краевого условия
    void third_boundary(double t)
    {
        std::fstream input;
        input.open("thirdboundary.txt");
        std::vector<int> BoundaryBorder;
        int variable;
        while (input >> variable)
        {
            BoundaryBorder.push_back(variable);
        }
        input.close();
        std::vector<double> LocUbeta = std::vector<double>(3);
        std::vector<std::vector<double>> LocBondA = std::vector<std::vector<double>>(3);
        for (int i = 0; i < 3; i++)
        {
            LocBondA[i].resize(3);
        }
        LocBondA[0][0] = 2. / 15.;
        LocBondA[0][1] = 1. / 15.;
        LocBondA[0][2] = -1. / 30.;
        LocBondA[1][0] = 1. / 15.;
        LocBondA[1][1] = 8. / 15.;
        LocBondA[1][2] = 1. / 15.;
        LocBondA[2][0] = -1. / 30.;
        LocBondA[2][1] = 1. / 15.;
        LocBondA[2][2] = 2. / 15.;
        for (int Border : BoundaryBorder)
        {
            switch (Border)
            {
            case 0:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1);
                    double bet = beta();
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _glob_b[_el_bonds[0][num][0]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][3]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][6]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 1:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i;
                    double bet = beta();
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _glob_b[_el_bonds[0][num][0]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][1]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][2]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 2:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1) + _nx - 2;
                    double bet = beta();
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _glob_b[_el_bonds[0][num][2]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][5]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][8]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 3:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i + (_nx - 1) * (_ny - 2);
                    double bet = beta();
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _glob_b[_el_bonds[0][num][6]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][7]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _glob_b[_el_bonds[0][num][8]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            }
        }
    }

    //Вставляет в глобальную матрицу локальную матрицу ребра при учёте третьего краевого условия
    void insert_loc_bond_A(int num, double h, int Border, std::vector<std::vector<double>> LocBondA)
    {
        switch (Border)
        {
        case 0:
        {
            double bet = beta();
            _di[_el_bonds[0][num][0]] += h * bet * LocBondA[0][0];
            _di[_el_bonds[0][num][3]] += h * bet * LocBondA[1][1];
            _di[_el_bonds[0][num][6]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_el_bonds[0][num][i * 3], _el_bonds[0][num][3 * j], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 1:
        {
            double bet = beta();
            _di[_el_bonds[0][num][0]] += h * bet * LocBondA[0][0];
            _di[_el_bonds[0][num][1]] += h * bet * LocBondA[1][1];
            _di[_el_bonds[0][num][2]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_el_bonds[0][num][i], _el_bonds[0][num][j], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 2:
        {
            double bet = beta();
            _di[_el_bonds[0][num][2]] += h * bet * LocBondA[0][0];
            _di[_el_bonds[0][num][5]] += h * bet * LocBondA[1][1];
            _di[_el_bonds[0][num][8]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_el_bonds[0][num][3 * i + 2], _el_bonds[0][num][3 * j + 2], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 3:
        {
            double bet = beta();
            _di[_el_bonds[0][num][6]] += h * bet * LocBondA[0][0];
            _di[_el_bonds[0][num][7]] += h * bet * LocBondA[1][1];
            _di[_el_bonds[0][num][8]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_el_bonds[0][num][6 + i], _el_bonds[0][num][6 + j], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        }
    }

    //Учёт первого краевого условия (Зануление строки)
    void first_boundary(double t)
    {
        std::fstream input;
        input.open("firstboundary.txt");
        std::vector<int> BoundaryBorder;
        int variable;
        while (input >> variable)
        {
            BoundaryBorder.push_back(variable);
        }
        input.close();

        for (int Border : BoundaryBorder)
        {
            switch (Border)
            {
            case 0:
            {
                int num;
                double C = 1e+30;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1);
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.0;
                    double Ugi;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_grid[_el_bonds[1][num][0]][0], _grid[_el_bonds[1][num][0]][1] + h * j, Border, t);
                        nullify_str(_el_bonds[0][num][3 * j], Ugi);
                        _di[_el_bonds[0][num][3 * j]] = 1.;
                        _glob_b[_el_bonds[0][num][3 * j]] = Ugi;
                    }
                }
                break;
            }
            case 1:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i;
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_grid[_el_bonds[1][num][0]][0] + h * j, _grid[_el_bonds[1][num][0]][1], Border, t);
                        nullify_str(_el_bonds[0][num][j], Ugi);
                        _di[_el_bonds[0][num][j]] = 1.;
                        _glob_b[_el_bonds[0][num][j]] = _di[_el_bonds[0][num][j]] * Ugi;
                    }
                }
                break;
            }
            case 2:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1) + _nx - 2;
                    double h = (_grid[_el_bonds[1][num][2]][1] - _grid[_el_bonds[1][num][0]][1]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_grid[_el_bonds[1][num][1]][0], _grid[_el_bonds[1][num][1]][1] + h * j, Border, t);
                        nullify_str(_el_bonds[0][num][3 * j + 2], Ugi);
                        _di[_el_bonds[0][num][3 * j + 2]] = 1.;
                        _glob_b[_el_bonds[0][num][3 * j + 2]] = _di[_el_bonds[0][num][3 * j + 2]] * Ugi;
                    }
                }
                break;
            }
            case 3:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i + (_nx - 1) * (_ny - 2);
                    double h = (_grid[_el_bonds[1][num][1]][0] - _grid[_el_bonds[1][num][0]][0]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_grid[_el_bonds[1][num][2]][0] + h * j, _grid[_el_bonds[1][num][2]][1], Border, t);
                        nullify_str(_el_bonds[0][num][6 + j], Ugi);
                        _di[_el_bonds[0][num][6 + j]] = 1.;
                        _glob_b[_el_bonds[0][num][6 + j]] = _di[_el_bonds[0][num][6 + j]] * Ugi;
                    }
                }
                break;
            }
            }
        }
        /*print_(_di, _au, _al);
        for (double x : _GlobB)
        {
           std::cout << x << std::endl;
        }*/
    }

    //Обнуление внедиагональных элементов в строке при учёте первого краевого условия
    void nullify_str(int Node, double Ugi)
    {
        int CurNode;
        for (int i = _ig[Node]; i < _ig[Node + 1]; i++)
        {
            _al[i] = 0.;
        }

        for (int i = Node; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            CurNode = _jg[_ig[i]];
            int k = 0;
            while ((CurNode <= Node) && (k < (_ig[i + 1] - _ig[i])))
            {
                if (CurNode == Node)
                {
                    _au[_ig[i] + k] = 0.;
                }
                k++;
                if (k < (_ig[i + 1] - _ig[i]))
                    CurNode = _jg[_ig[i] + k];
            }
        }
    }

    //Факторизация матрицы методом LUsq для решения СЛАУ
    void LUsq(std::vector<double> &_dif, std::vector<double> &ggl_f, std::vector<double> &ggu_f)
    {
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            double sumdi = 0.0;

            int i0 = _ig[i];
            int i1 = _ig[i + 1];

            for (int k = i0; k < i1; k++)
            {
                int j = _jg[k];
                int j0 = _ig[j];

                int j1 = _ig[j + 1];

                int ik = i0;
                int kj = j0;

                double suml = 0.0;
                double sumu = 0.0;

                while (ik < k)
                {

                    if (_jg[ik] == _jg[kj])
                    {

                        suml += ggl_f[ik] * ggu_f[kj];
                        sumu += ggu_f[ik] * ggl_f[kj];
                        ik++;
                        kj++;
                    }

                    else
                        _jg[ik] > _jg[kj] ? kj++ : ik++;
                }

                ggl_f[k] = (ggl_f[k] - suml) / _dif[j];
                ggu_f[k] = (ggu_f[k] - sumu) / _dif[j];
                sumdi += ggl_f[k] * ggu_f[k];
            }

            _dif[i] = sqrt(_dif[i] - sumdi);
        }
        // print_(_dif, ggu_f, ggl_f);
    }

    //Умножение матрицы на вектор
    std::vector<double> mult(std::vector<double> &v)
    {
        std::vector<double> res(v.size());
        for (int i = 0; i < v.size(); i++)
        {
            res[i] = _di[i] * v[i];
            for (int j = _ig[i]; j < _ig[i + 1]; j++)
            {
                res[i] += _al[j] * v[_jg[j]];
                res[_jg[j]] += _au[j] * v[i];
            }
        }
        return res;
    }

    //Решение СЛАУ с помощью локально-оптимальной схемы с предобуславливанием
    void LoS_precond(int j)
    {
        std::vector<double> _dif = _di;
        std::vector<double> _auf = _au;
        std::vector<double> _alf = _al;
        LUsq(_dif, _alf, _auf);
        int k = 0;
        std::vector<double> buf = mult(_q[j]);
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            buf[i] = _glob_b[i] - buf[i];
        }
        std::vector<double> r = LU_direct(buf, _dif, _alf);
        double error = scalar_prod(r, r);
        double error1 = error + 1;
        std::vector<double> z = LU_reverse(r, _dif, _auf);
        buf = mult(z);
        std::vector<double> p = LU_direct(buf, _dif, _alf);
        ;
        while (error > 1e-15 && k < 1000 && abs((error - error1)) >= 1e-16)
        {
            double pp = scalar_prod(p, p);
            double pr = scalar_prod(p, r);
            double alpha = pr / pp;
            error1 = error;
            error -= alpha * alpha * pp;
            for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
            {
                _q[j][i] += alpha * z[i];
                r[i] -= alpha * p[i];
            }
            std::vector<double> Ur = LU_reverse(r, _dif, _auf);
            buf = mult(Ur);
            buf = LU_direct(buf, _dif, _alf);
            double betta = -(scalar_prod(p, buf) / pp);
            for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
            {
                z[i] = Ur[i] + betta * z[i];
                p[i] = buf[i] + betta * p[i];
            }
            k++;
        }
        std::cout << "k:" << k << "\nerror: " << error << "\n";
        /*for (int i = 0; i < _Q.size(); i++)
        {
           std::cout << _Q[i] << std::endl;
        }*/
    };

    //Прямой ход для ЛОС
    std::vector<double> LU_direct(const std::vector<double> &b, std::vector<double> &_dif, std::vector<double> &_alf)
    {
        std::vector<double> res = b;

        for (size_t i = 0; i < res.size(); i++)
        {
            double sum = 0.0;
            for (size_t j = _ig[i]; j < _ig[i + 1]; j++)
                sum += _alf[j] * res[_jg[j]];
            res[i] -= sum;
            res[i] /= _dif[i];
        }
        return res;
    }

    //Обратный ход для ЛОС
    std::vector<double> LU_reverse(const std::vector<double> &b, std::vector<double> &_dif, std::vector<double> &ggu_f)
    {
        std::vector<double> res = b;

        for (int i = (2 * _nx - 1) * (2 * _ny - 1) - 1; i >= 0; i--)
        {
            res[i] /= _dif[i];
            for (size_t j = _ig[i]; j < _ig[i + 1]; j++)
                res[_jg[j]] -= ggu_f[j] * res[i];
        }
        return res;
    }

    //Скалярное произведение двух векторов
    double scalar_prod(std::vector<double> x, std::vector<double> y)
    {
        double res = 0.0;
        if (x.size() == y.size())
        {
            for (int i = 0; i < y.size(); i++)
            {
                res += x[i] * y[i];
            }
            return res;
        }
        else
        {
            std::cout << "Error!";
            return res;
        }
    }

    //Вывод матрицы для возможных проверок
    void print_mx(std::vector<double> di, std::vector<double> ggu, std::vector<double> ggl)
    {
        std::vector<std::vector<double>> mat;
        mat.resize((2 * _nx - 1) * (2 * _ny - 1));
        for (int i = 0; i < mat.size(); i++)
        {
            mat[i].resize((2 * _nx - 1) * (2 * _ny - 1));
        }

        for (int i = 0; i < mat.size(); i++)
        {
            mat[i][i] = di[i];
            for (int j = _ig[i]; j < _ig[i + 1]; j++)
            {
                mat[i][_jg[j]] = ggl[j];
                mat[_jg[j]][i] = ggu[j];
            }
        }

        for (int i = 0; i < mat.size(); i++)
        {
            for (int j = 0; j < mat[i].size(); j++)
                std::cout << std::setw(12) << mat[i][j] << " ";
            std::cout << std::endl;
        }
    }

    //Вывод правой части СЛАУ для проверки
    void print_b()
    {
        int n = _glob_b.size();
        std::cout << std::endl
                  << "RhsVec" << std::endl;
        for (int i = 0; i < n; i++)
        {
            std::cout << _glob_b[i] << std::endl;
        }
    }

    //Выводит решение СЛАУ, значение искомой функции в узлах и их разницу
    void print_error()
    {
        double hx, hy, x, y;
        int size = _el_bonds[0].size();
        std::vector<std::vector<double>> fictive_grid = std::vector<std::vector<double>>((2 * _nx - 1) * (2 * _ny - 1));
        for (int i = 0; i < size; i++)
        {
            x = _grid[_el_bonds[1][i][0]][0];
            y = _grid[_el_bonds[1][i][0]][1];
            hx = (_grid[_el_bonds[1][i][1]][0] - _grid[_el_bonds[1][i][0]][0]) / 2.;
            hy = (_grid[_el_bonds[1][i][2]][1] - _grid[_el_bonds[1][i][0]][1]) / 2.;
            for (int j = 0; j < 9; j++)
            {
                fictive_grid[_el_bonds[0][i][j]].resize(2);
                fictive_grid[_el_bonds[0][i][j]][0] = x + hx * (j % 3);
                fictive_grid[_el_bonds[0][i][j]][1] = y + hy * (j / 3);
            }
        }

        size = fictive_grid.size();

        std::cout << std::setw(36) << "(x;y;t)"
                  << " | " << std::setw(12) << "q(i)"
                  << " | " << std::setw(12)
                  << "u(i)"
                  << " | " << std::setw(14)
                  << "|q(i) - u(i)|"
                  << " | " << std::endl;
        std::cout << "_____________________________________|______________|______________|________________|" << std::endl;
        int nt = _time.size();
        double tolerance = 0.;
        double real_u = 0.;
        for (int j = 0; j < nt; j++)
        {
            for (int i = 0; i < size; i++)
            {
                std::cout << "(" << std::setw(8) << fictive_grid[i][0] << ";"
                          << std::setw(12) << fictive_grid[i][1] << ";"
                          << std::setw(12) << _time[j] << ")"
                                                          " | "
                          << std::setw(12) << _q[j][i] << " | " << std::setw(12)
                          << U(fictive_grid[i][0], fictive_grid[i][1], _time[j]) << " | " << std::setw(14) << abs(_q[j][i] - U(fictive_grid[i][0], fictive_grid[i][1], _time[j])) << " | " << std::endl;
                tolerance += abs(_q[j][i] - U(fictive_grid[i][0], fictive_grid[i][1], _time[j])) * abs(_q[j][i] - U(fictive_grid[i][0], fictive_grid[i][1], _time[j]));
                real_u += U(fictive_grid[i][0], fictive_grid[i][1], _time[j]) * U(fictive_grid[i][0], fictive_grid[i][1], _time[j]);
            }
        }
        std::cout << "Relative tolerance = " << sqrt(tolerance / real_u) << std::endl
                  << "Tolerance = " << sqrt(tolerance) << std::endl;
        /*double k = 0;
        for (int i = 0; i < 2 * _ny - 1; i++)
        {
           for (int l = 0; l < 2 * _nx - 1; l++)
           {
              double x = _grid[_ElBonds[1][0][0]][0] + hx * l;
              double y = _grid[_ElBonds[1][0][0]][1] + hy * i;
              std::cout << "(" << std::setw(7) << x << ";"
                 << std::setw(7) << y << ")"
                 " | " << std::setw(12) << _Q[j][k] << " | " << std::setw(12)
                 << U(x, y,0) << " | " << std::setw(12) <<
                 _Q[j][k] - U(x, y, 0) << " | " << std::endl;
              k++;
           }
        }*/
    }

    //Решаем гиперболическую задачу
    void solve_hyperbolic()
    {
        making_time_grid();
        making_initial_conds();
        three_layer_global_bild(_time[0], _time[1], _time[2]);
        second_boundary(_time[2]);
        third_boundary(_time[2]);
        first_boundary(_time[2]);
        init_hyperbolic_approx(2);
        LoS_precond(2);
        print_error();
    }
};

int main()
{
    slae SLAE;
    int method;
    SLAE.making_grid(); //Загружаем область решения и разбиваем её на конечные элементы
    SLAE.solve_hyperbolic();
    return 0;
}

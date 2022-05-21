#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

//ТУТ ПОЯСНИТЕЛЬНАЯ ЗАПИСКА//
/*В файле input.txt задаётся расчётная область: граничные точки по х, граничные точки по у, количество узлов по х, количество узлов по у
  В файлах firstboundary.txt, secondboundary.txt, thirdboundary.txt находятся номера границ, на которых заданы соответствующие краевые условия
   Для учёта краевых условий пронумеровали от 0 до 3 границы расчётной области против часовой стрелки начиная с левой границы
   Функции и коэффициенты для краевых условий задаются внутри программы в соответствующих функциях для нужных границ
   Формат хранения глобальной матрицы - разреженный*/

class slae
{
private:
    int _nx;              // количество узлов по оси х
    int _ny;              // количество узлов по оси у
    std::vector<int> _ig; // Профиль глобальной матрицы
    std::vector<int> _jg;
    std::vector<double> _al;                             //Нижний треугольник матрицы
    std::vector<double> _au;                             //Верхний треугольник матрицы
    std::vector<double> _di;                             //Диагональные элементы матрицы
    std::vector<std::vector<double>> _LocMass;           //Локальная матрица масс
    std::vector<std::vector<std::vector<double>>> _LocG; // Локальная матрицы жёсткости
    std::vector<std::vector<double>> _Grid;              //Пронумерованные точки в узлах сетки
    std::vector<double> _time;                           //Сетка по времени
    std::vector<std::vector<std::vector<int>>> _ElBonds; // Связи узлов и элементов
    std::vector<std::vector<int>> _NodeBonds;
    std::vector<double> _GlobB;          //Вектор правой части
    std::vector<std::vector<double>> _Q; //Вектор решения СЛАУ (коэффициенты разложения искомой функции по биквадратичному базису)
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
        this->_LocMass = std::vector<std::vector<double>>();
        this->_LocG = std::vector<std::vector<std::vector<double>>>();
        this->_Grid = std::vector<std::vector<double>>();
        this->_ElBonds = std::vector<std::vector<std::vector<int>>>();
        this->_NodeBonds = std::vector<std::vector<int>>();
        this->_GlobB = std::vector<double>();
        this->_Q = std::vector<std::vector<double>>();
        this->_time = std::vector<double>();
    }
    /*void LoadData()
    {
       std::fstream Grid;
       Grid.open("Grid.txt");
       Grid >> _NodeNum;
       _Grid.resize(_NodeNum);
       for (int i = 0; i < _NodeNum; i++)
       {
          _Grid[i].resize(2);
          Grid >> _Grid[i][0] >> _Grid[i][1];
       }
       Grid.close();
    }*/

    //Вовзращает значение функции правой части уравнение в точке
    double f(double x, double y, double t)
    {
        return 1.;
    }

    //Возвращает значение функции коэффициента диффузии в точке
    double lambda(double x, double y)
    {
        return 1.;
    }

    //Возвращает значение коэффициента гамма
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
            return t;
            break;
        }
        case 1:
        {
            return t;
            break;
        }
        case 2:
        {
            return t;
            break;
        }
        case 3:
        {
            return t;
            break;
        }
        }
    }

    //Возвращает значение искомой функции в точке для проверки точности решения
    double U(double x, double y, double t)
    {
        return t;
    }

    //Начальное значение
    double U_0(double x, double y)
    {
        return 1.;
    }

    //Начальная производная
    double dU_0(double x, double y)
    {
        return 1.0;
    }

    double U_1(double x, double y)
    {
        return x * x + 1.010025;
    }

    //Создаёт сетку и список связи элементов и узлов
    void making_grid()
    {
        std::fstream input;
        input.open("input.txt");
        double xmin, xmax, ymin, ymax;
        input >> xmin;
        input >> xmax;
        input >> ymin;
        input >> ymax;
        input >> _nx;
        input >> _ny;
        double hx = (xmax - xmin) / (_nx - 1);
        double hy = (ymax - ymin) / (_ny - 1);

        _Grid.resize(_nx * _ny);
        int k = 0;
        for (int i = 0; i < _ny; i++)
            for (int j = 0; j < _nx; j++)
            {
                _Grid[k].resize(2);
                _Grid[k][0] = xmin + j * hx;
                _Grid[k][1] = ymin + i * hy;
                k++;
            }
        input.close();
        _ElBonds.resize(2);
        _ElBonds[0].resize((_nx - 1) * (_ny - 1));
        _ElBonds[1].resize((_nx - 1) * (_ny - 1));
        for (int i = 0; i < (_ny - 1); i++)
        {
            for (int j = 0; j < (_nx - 1); j++)
            {
                _ElBonds[0][i * (_nx - 1) + j].resize(9);
                _ElBonds[0][i * (_nx - 1) + j][0] = 2 * i * (2 * _nx - 1) + j * 2;
                _ElBonds[0][i * (_nx - 1) + j][1] = 2 * i * (2 * _nx - 1) + j * 2 + 1;
                _ElBonds[0][i * (_nx - 1) + j][2] = 2 * i * (2 * _nx - 1) + j * 2 + 2;
                _ElBonds[0][i * (_nx - 1) + j][3] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx - 1;
                _ElBonds[0][i * (_nx - 1) + j][4] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx;
                _ElBonds[0][i * (_nx - 1) + j][5] = 2 * i * (2 * _nx - 1) + j * 2 + 2 * _nx + 1;
                _ElBonds[0][i * (_nx - 1) + j][6] = 2 * (i + 1) * (2 * _nx - 1) + j * 2;
                _ElBonds[0][i * (_nx - 1) + j][7] = 2 * (i + 1) * (2 * _nx - 1) + j * 2 + 1;
                _ElBonds[0][i * (_nx - 1) + j][8] = 2 * (i + 1) * (2 * _nx - 1) + j * 2 + 2;
                _ElBonds[1][i * (_nx - 1) + j].resize(4);
                _ElBonds[1][i * (_nx - 1) + j][0] = i * _nx + j;
                _ElBonds[1][i * (_nx - 1) + j][1] = i * _nx + j + 1;
                _ElBonds[1][i * (_nx - 1) + j][2] = (i + 1) * _nx + j;
                _ElBonds[1][i * (_nx - 1) + j][3] = (i + 1) * _nx + j + 1;
            }
        }

        _NodeBonds.resize((2 * _nx - 1) * (2 * _ny - 1));
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            for (int j = 0; j < 9; j++)
            {
                for (int k = 0; k < 9; k++)
                {
                    if (_ElBonds[0][i][j] > _ElBonds[0][i][k])
                    {
                        _NodeBonds[_ElBonds[0][i][j]].push_back(_ElBonds[0][i][k]);
                    }
                }
            }
        }

        for (int i = 0; i < ((2 * _nx - 1) * (2 * _ny - 1)); i++)
        {
            sort(_NodeBonds[i].begin(), _NodeBonds[i].end());
            auto last = std::unique(_NodeBonds[i].begin(), _NodeBonds[i].end());
            _NodeBonds[i].erase(last, _NodeBonds[i].end());
        }
    }

    //Создаём сетку по времени
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
        _Q.resize(nt);
    }

    void making_initial_conds()
    {
        double xmin, ymin, hx, hy;
        _Q[0].resize((2 * _nx - 1) * (2 * _ny - 1));
        _Q[1].resize((2 * _nx - 1) * (2 * _ny - 1));
        xmin = _Grid[0][0];
        ymin = _Grid[0][1];
        hx = (_Grid[1][0] - _Grid[0][0]) / 2.;
        hy = (_Grid[_nx][1] - _Grid[0][1]) / 2.;
        int k = 0;
        for (int i = 0; i < 2 * _ny - 1; i++)
        {
            for (int j = 0; j < 2 * _nx - 1; j++)
            {
                _Q[0][k] = U_0(xmin + hx * j, ymin + hy * i);
                _Q[1][k] = _Q[0][k] + dU_0(xmin + hx * j, ymin + hy * i) * (_time[1] - _time[0]);
                //_Q[1][k] = U_1(xmin + hx * j, ymin + hy * i);
                k++;
            }
        }
    }

    void init_elliptic_approx()
    {
        _Q.resize(1);
        _Q[0].resize((2 * _nx - 1) * (2 * _ny - 1));
    }

    void init_hyperbolic_approx(int j)
    {
        _Q[j].resize((2 * _nx - 1) * (2 * _ny - 1));
        std::fill(_Q[j].begin(), _Q[j].end(), 1.);
    }

    //Создаёт профиль глобальной матрицы
    void set_profile()
    {
        _ig.resize((2 * _nx - 1) * (2 * _ny - 1) + 1);
        _di.resize((2 * _nx - 1) * (2 * _ny - 1));
        _ig[0] = 0;
        _ig[1] = 0;
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            int k = 0;
            for (int j = 0; j < _NodeBonds[i].size(); j++)
            {
                k++;
                _jg.push_back(_NodeBonds[i][j]);
            }
            _ig[i + 1] = _ig[i] + k;
        }
    }

    //Производит сборку локальных матриц
    void local_build()
    {
        std::vector<std::vector<std::vector<double>>> LocMassPsi = std::vector<std::vector<std::vector<double>>>(3);
        std::vector<std::vector<std::vector<double>>> LocGPsi = std::vector<std::vector<std::vector<double>>>(2);
        double hx, hy;
        hx = (_Grid[1][0] - _Grid[0][0]);
        hy = (_Grid[_nx][1] - _Grid[0][1]);
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

        for (int i = 0; i < 2; i++)
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

        _LocMass.resize(9);
        for (int i = 0; i < 9; i++)
        {
            _LocMass[i].resize(9);
            for (int j = 0; j < 9; j++)
            {
                _LocMass[i][j] = hx * hy * LocMassPsi[2][i / 3][j / 3] * LocMassPsi[2][i % 3][j % 3];
            }
        }
        _LocG.resize(4);
        for (int k = 0; k < 4; k++)
        {
            _LocG[k].resize(9);
        }
        for (int i = 0; i < 9; i++)
        {
            _LocG[0][i].resize(9);
            _LocG[1][i].resize(9);
            _LocG[2][i].resize(9);
            _LocG[3][i].resize(9);
            for (int j = 0; j < 9; j++)
            {
                _LocG[0][i][j] = hy * LocGPsi[0][i % 3][j % 3] * LocMassPsi[0][i / 3][j / 3] / hx + hx * LocGPsi[0][i / 3][j / 3] * LocMassPsi[0][i % 3][j % 3] / hy;
                _LocG[1][i][j] = hy * LocGPsi[0][i % 3][j % 3] * LocMassPsi[1][i / 3][j / 3] / hx + hx * LocGPsi[1][i / 3][j / 3] * LocMassPsi[0][i % 3][j % 3] / hy;
                _LocG[2][i][j] = hy * LocGPsi[1][i % 3][j % 3] * LocMassPsi[0][i / 3][j / 3] / hx + hx * LocGPsi[0][i / 3][j / 3] * LocMassPsi[1][i % 3][j % 3] / hy;
                _LocG[3][i][j] = hy * LocGPsi[1][i % 3][j % 3] * LocMassPsi[1][i / 3][j / 3] / hx + hx * LocGPsi[1][i / 3][j / 3] * LocMassPsi[1][i % 3][j % 3] / hy;
            }
        }
    }

    //Производит сборку глобальной матрицы по элементам
    void elliptic_global_build(double t)
    {
        double gamma = 0.;
        std::vector<double> LocB = std::vector<double>(9);
        std::vector<std::vector<double>> LocG;
        _GlobB.resize((2 * _nx - 1) * (2 * _ny - 1));
        _al.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        _au.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            get_local_gamma(i, gamma);
            get_local_f(LocB, i, t);
            mass_mx_mult_vec(LocB);
            for (int k = 0; k < 9; k++)
            {
                _GlobB[_ElBonds[0][i][k]] += LocB[k];
            }
            get_loc_G(i, LocG);
            for (int k = 0; k < 9; k++)
            {
                _di[_ElBonds[0][i][k]] += gamma * _LocMass[k][k] + LocG[k][k];
            }

            int Index;
            for (int k = 1; k < 9; k++)
            {
                for (int j = 0; j < k; j++)
                {
                    get_index(_ElBonds[0][i][k], _ElBonds[0][i][j], Index);
                    _al[Index] += gamma * _LocMass[k][j] + LocG[k][j];
                    _au[Index] += gamma * _LocMass[j][k] + LocG[j][k];
                }
            }
        }
    }

    void three_layer_global_build(double t0, double t1, double t2)
    {
        double gamma = 0.;
        double sigma = 0.;
        double khi = 0.;

        std::vector<double> loc_b = std::vector<double>(9);
        std::vector<double> q_0 = std::vector<double>(9);
        std::vector<double> q_1 = std::vector<double>(9);
        std::vector<double> q_1_G = std::vector<double>(9);
        std::vector<std::vector<double>> LocG;

        double dt = t2 - t0;
        double dt0 = t2 - t1;
        double dt1 = t1 - t0;

        _GlobB.resize((2 * _nx - 1) * (2 * _ny - 1));
        _al.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        _au.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            get_local_gamma(i, gamma);
            get_local_sigma(i, sigma);
            get_local_khi(i, khi);
            get_local_f(loc_b, i, t1);

            for (int k = 0; k < 9; k++)
            {
                q_0[k] = _Q[0][_ElBonds[0][i][k]];
                q_1[k] = _Q[1][_ElBonds[0][i][k]];
            }

            q_1_G = q_1;
            mass_mx_mult_vec(loc_b);
            mass_mx_mult_vec(q_0);
            mass_mx_mult_vec(q_1);
            get_loc_G(i, LocG);
            stiffness_mx_mult_vec(q_1_G, LocG);

            for (int k = 0; k < 9; k++)
            {
                _GlobB[_ElBonds[0][i][k]] += loc_b[k] - khi * (2. * q_0[k] / (dt * dt1) - 2. * q_1[k] / (dt1 * dt0)) -
                                             sigma * ((dt0 - dt1) * q_1[k] / (dt1 * dt0) - dt0 * q_0[k] / (dt * dt1)) -
                                             gamma * q_1[k] - q_1_G[k];
            }

            for (int k = 0; k < 9; k++)
            {
                _di[_ElBonds[0][i][k]] += 2. * khi * _LocMass[k][k] / (dt * dt0) + sigma * dt1 * _LocMass[k][k] / (dt0 * dt);
            }

            int Index;
            for (int k = 1; k < 9; k++)
            {
                for (int j = 0; j < k; j++)
                {
                    get_index(_ElBonds[0][i][k], _ElBonds[0][i][j], Index);
                    _al[Index] += 2. * khi * _LocMass[k][j] / (dt * dt0) + sigma * dt1 * _LocMass[k][j] / (dt0 * dt);
                    _au[Index] += 2. * khi * _LocMass[k][j] / (dt * dt0) + sigma * dt1 * _LocMass[k][j] / (dt0 * dt);
                }
            }
        }
    }

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
                LocG[i][j] = LocPhiLam[0] * _LocG[0][i][j] + LocPhiLam[1] * _LocG[1][i][j] + LocPhiLam[2] * _LocG[2][i][j] + LocPhiLam[3] * _LocG[3][i][j];
            }
        }
    }

    //Возвращает коэффициенты разложения коэффициента диффузии по билинейным базисным функциям
    void get_loc_lam(int num, std::vector<double> &LocPhiLam)
    {
        for (int i = 0; i < 4; i++)
        {
            LocPhiLam[i] = lambda(_Grid[_ElBonds[1][num][i]][0], _Grid[_ElBonds[1][num][i]][1]);
        }
        int k = 0;
    }

    //Возвращают усреднённые коэффициенты при матрице масс для элемента
    void get_local_gamma(int num, double &gam)
    {
        double hx, hy;
        double g = 0.0;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            g += gamma(_Grid[_ElBonds[1][num][0]][0] + hx * (i % 3), _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3));
        }
        gam = g / 9.;
    }

    void get_local_sigma(int num, double &sig)
    {
        double hx, hy;
        double s = 0.0;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            s += sigma(_Grid[_ElBonds[1][num][0]][0] + hx * (i % 3), _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3));
        }
        sig = s / 9.;
    }

    void get_local_khi(int num, double &kh)
    {
        double hx, hy;
        double k = 0.0;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            k += khi(_Grid[_ElBonds[1][num][0]][0] + hx * (i % 3), _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3));
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
                sum += f[i] * _LocMass[i][j];
            }
            res[i] = sum;
        }
        f = res;
    }

    //Получает вектор локальных значений правой части на конечном элементе
    void get_local_f(std::vector<double> &LocB, int num, double t)
    {
        double hx, hy;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        std::vector<double> res = std::vector<double>(9);
        for (int i = 0; i < 9; i++)
        {
            double x = _Grid[_ElBonds[1][num][0]][0] + hx * (i % 3);
            double y = _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3);
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    get_local_theta(num, h / 2.0, Border, LocTheta, t);
                    _GlobB[_ElBonds[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][3]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][6]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _GlobB[_ElBonds[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][1]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][2]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _GlobB[_ElBonds[0][num][2]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][5]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    get_local_theta(num, h / 2., Border, LocTheta, t);
                    _GlobB[_ElBonds[0][num][6]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][7]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
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
                LocTheta[i] = theta(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * i, Border, t);
                break;
            }
            case 1:
            {
                LocTheta[i] = theta(_Grid[_ElBonds[1][num][0]][0] + h * i, _Grid[_ElBonds[1][num][0]][1], Border, t);
                break;
            }
            case 2:
            {
                LocTheta[i] = theta(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * i, Border, t);
                break;
            }
            case 3:
            {
                LocTheta[i] = theta(_Grid[_ElBonds[1][num][2]][0] + h * i, _Grid[_ElBonds[1][num][2]][1], Border, t);
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
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * i, Border, t);
                break;
            }
            case 1:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][0]][0] + h * i, _Grid[_ElBonds[1][num][0]][1], Border, t);
                break;
            }
            case 2:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * i, Border, t);
                break;
            }
            case 3:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][2]][0] + h * i, _Grid[_ElBonds[1][num][2]][1], Border, t);
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _GlobB[_ElBonds[0][num][0]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][3]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][6]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _GlobB[_ElBonds[0][num][0]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][1]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][2]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _GlobB[_ElBonds[0][num][2]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][5]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    insert_loc_bond_A(num, h, Border, LocBondA);
                    get_local_Ubeta(num, h / 2., Border, LocUbeta, t);
                    _GlobB[_ElBonds[0][num][6]] += h * bet * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][7]] += h * bet * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * bet * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
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
            _di[_ElBonds[0][num][0]] += h * bet * LocBondA[0][0];
            _di[_ElBonds[0][num][3]] += h * bet * LocBondA[1][1];
            _di[_ElBonds[0][num][6]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_ElBonds[0][num][i * 3], _ElBonds[0][num][3 * j], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 1:
        {
            double bet = beta();
            _di[_ElBonds[0][num][0]] += h * bet * LocBondA[0][0];
            _di[_ElBonds[0][num][1]] += h * bet * LocBondA[1][1];
            _di[_ElBonds[0][num][2]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_ElBonds[0][num][i], _ElBonds[0][num][j], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 2:
        {
            double bet = beta();
            _di[_ElBonds[0][num][2]] += h * bet * LocBondA[0][0];
            _di[_ElBonds[0][num][5]] += h * bet * LocBondA[1][1];
            _di[_ElBonds[0][num][8]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_ElBonds[0][num][3 * i + 2], _ElBonds[0][num][3 * j + 2], Index);
                    _al[Index] += h * bet * LocBondA[i][j];
                    _au[Index] += h * bet * LocBondA[i][j];
                }
            }
            break;
        }
        case 3:
        {
            double bet = beta();
            _di[_ElBonds[0][num][6]] += h * bet * LocBondA[0][0];
            _di[_ElBonds[0][num][7]] += h * bet * LocBondA[1][1];
            _di[_ElBonds[0][num][8]] += h * bet * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    get_index(_ElBonds[0][num][6 + i], _ElBonds[0][num][6 + j], Index);
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.0;
                    double Ugi;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * j, Border, t);
                        nullify_str(_ElBonds[0][num][3 * j], Ugi);
                        _di[_ElBonds[0][num][3 * j]] = 1.;
                        _GlobB[_ElBonds[0][num][3 * j]] = Ugi;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][0]][0] + h * j, _Grid[_ElBonds[1][num][0]][1], Border, t);
                        nullify_str(_ElBonds[0][num][j], Ugi);
                        _di[_ElBonds[0][num][j]] = 1.;
                        _GlobB[_ElBonds[0][num][j]] = _di[_ElBonds[0][num][j]] * Ugi;
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
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * j, Border, t);
                        nullify_str(_ElBonds[0][num][3 * j + 2], Ugi);
                        _di[_ElBonds[0][num][3 * j + 2]] = 1.;
                        _GlobB[_ElBonds[0][num][3 * j + 2]] = _di[_ElBonds[0][num][3 * j + 2]] * Ugi;
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
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][2]][0] + h * j, _Grid[_ElBonds[1][num][2]][1], Border, t);
                        nullify_str(_ElBonds[0][num][6 + j], Ugi);
                        _di[_ElBonds[0][num][6 + j]] = 1.;
                        _GlobB[_ElBonds[0][num][6 + j]] = _di[_ElBonds[0][num][6 + j]] * Ugi;
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
        std::vector<double> buf = mult(_Q[j]);
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            buf[i] = _GlobB[i] - buf[i];
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
                _Q[j][i] += alpha * z[i];
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
    double scalar_prod(std::vector<double> &x, std::vector<double> &y)
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

    void print_b()
    {
        int n = _GlobB.size();
        std::cout << std::endl
                  << "RhsVec" << std::endl;
        for (int i = 0; i < n; i++)
        {
            std::cout << _GlobB[i] << std::endl;
        }
    }

    //Выводит решение СЛАУ, значение искомой функции в узлах и их разницу
    void print_error(int j)
    {
        double hx, hy;
        hx = (_Grid[_ElBonds[1][0][1]][0] - _Grid[_ElBonds[1][0][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][0][2]][1] - _Grid[_ElBonds[1][0][0]][1]) / 2.;
        std::cout << "t = " << _time[j] << std::endl;
        std::cout << std::setw(17) << "(x;y)"
                  << " | " << std::setw(12) << "q(i)"
                  << " | " << std::setw(12)
                  << "u(i)"
                  << " | " << std::setw(12)
                  << "q(i) - u(i)"
                  << " | " << std::endl;
        std::cout << "__________________|______________|______________|______________|" << std::endl;
        // for (int i = 0; i < (2*_nx - 1)*(2*_ny - 1); i++)
        //{
        //    double x = _Grid[_ElBonds[1][0][0]][0] + hx * (i % (2 * _nx - 1));
        //    double y = _Grid[_ElBonds[1][0][0]][1] + hy * (i / (2 * _ny - 1));
        //    std::cout << "(" << std::setw(7) << x << ";"
        //       << std::setw(7) << y << ")"
        //       " | " << std::setw(12) << _Q[j][i] << " | " << std::setw(12)
        //       << U(x, y, _time[j]) << " | " << std::setw(12) <<
        //          _Q[j][i] - U(x, y, _time[j]) << " | " << std::endl;
        // }

        double k = 0;
        for (int i = 0; i < 2 * _ny - 1; i++)
        {
            for (int l = 0; l < 2 * _nx - 1; l++)
            {
                double x = _Grid[_ElBonds[1][0][0]][0] + hx * l;
                double y = _Grid[_ElBonds[1][0][0]][1] + hy * i;
                std::cout << "(" << std::setw(7) << x << ";"
                          << std::setw(7) << y << ")"
                                                  " | "
                          << std::setw(12) << _Q[j][k] << " | " << std::setw(12)
                          << U(x, y, _time[j]) << " | " << std::setw(12) << _Q[j][k] - U(x, y, _time[j]) << " | " << std::endl;
                k++;
            }
        }
    }

    //Решаем гиперболическую задачу
    void solve_hyperbolic()
    {
        making_time_grid();
        making_initial_conds();
        three_layer_global_build(_time[0], _time[1], _time[2]);
        // print_mx(_di, _au, _al);
        // print_b();
        second_boundary(_time[2]);
        third_boundary(_time[2]);
        first_boundary(_time[2]);
        // print_mx(_di, _au, _al);
        // print_b();
        init_hyperbolic_approx(2);
        LoS_precond(2);
        three_layer_global_build(_time[1], _time[2], _time[3]);
        second_boundary(_time[3]);
        third_boundary(_time[3]);
        first_boundary(_time[3]);
        init_hyperbolic_approx(3);
        LoS_precond(3);
        print_error(2);
        print_error(3);
        /*int nt = _time.size();
        for (int j = 3; j < nt; j++)
        {
           four_layer_global_build(j);
           second_boundary(_time[j]);
           third_boundary(_time[j]);
           first_boundary(_time[j]);
           init_hyperbolic_approx(j);
           LoS_precond(j);
        }

        for (int j = 0; j < nt; j++)
        {
           print_error(j);
        }*/
    }

    //Решаем
};

int main()
{
    slae SLAE;
    int method;
    SLAE.making_grid(); //Загружаем область решения и разбиваем её на конечные элементы
    SLAE.set_profile(); //Инициализируем профиль глобальной матрицы
    SLAE.local_build(); //Собираем локальную матрицу
    SLAE.solve_hyperbolic();
}
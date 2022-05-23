#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cmath>

class matrix
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
    std::vector<std::vector<std::vector<int>>> _ElBonds; // Связи узлов и элементов
    std::vector<std::vector<int>> _NodeBonds;
    std::vector<double> _GlobB; //Вектор правой части
    std::vector<double> _Q;     //Вектор решения СЛАУ (коэффициенты разложения искомой функции по биквадратичному базису)
public:
    matrix()
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
        this->_Q = std::vector<double>();
    }
    matrix(int nx,
           int ny,
           std::vector<int> ig,
           std::vector<int> jg,
           std::vector<double> al,
           std::vector<double> au,
           std::vector<double> di,
           std::vector<std::vector<double>> LocMass,
           std::vector<std::vector<std::vector<double>>> LocG,
           std::vector<std::vector<double>> Grid,
           std::vector<std::vector<std::vector<int>>> ElBonds,
           std::vector<std::vector<int>> NodeBonds,
           std::vector<double> GlobB,
           std::vector<double> Q)
    {
        _nx = nx;
        _ny = ny;
        _ig = ig;
        _jg = jg;
        _al = al;
        _au = au;
        _di = di;
        _LocMass = LocMass;
        _LocG = LocG;
        _Grid = Grid;
        _ElBonds = ElBonds;
        _NodeBonds = NodeBonds;
        _GlobB = GlobB;
        _Q = Q;
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
    double F(double x, double y)
    {
        return 0.;
    }

    //Возвращает значение функции коэффициента диффузии в точке
    double Lambda(double x, double y)
    {
        return 1.;
    }

    //Возвращает значение коэффициента гамма
    double Gamma(double x, double y)
    {
        return 0.;
    }

    //Возвращает значение функции Тета во втором краевом условии
    double Theta(double x, double y, int Border)
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
    double Beta()
    {
        return 0.5;
    }

    //Возвращает значение функции Убета В третьем краевом условии
    double Ubeta(double x, double y, int Border)
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
    double Ug(double x, double y, int Border)
    {
        switch (Border)
        {
        case 0:
        {
            return 1. + y * y;
            break;
        }
        case 1:
        {
            return 1;
            break;
        }
        case 2:
        {
            return 2.;
            break;
        }
        case 3:
        {
            return 9.;
            break;
        }
        }
    }

    //Возвращает значение искомой функции в точке для проверки точности решения
    double U(double x, double y)
    {

        return 1.;
    }

    //Создаёт сетку и список связи элементов и узлов
    void makingGrid()
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

    //Создаёт профиль глобальной матрицы
    void SetProfile()
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
        int u = 0;
    }

    //Производит сборку локальных матриц
    void LocalBuild()
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
    void GlobalBuild()
    {
        double gamma = 0;
        std::vector<double> LocB = std::vector<double>(9);
        std::vector<std::vector<double>> LocG;
        _GlobB.resize((2 * _nx - 1) * (2 * _ny - 1));
        _al.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        _au.resize(_ig[(2 * _nx - 1) * (2 * _ny - 1)]);
        for (int i = 0; i < (_nx - 1) * (_ny - 1); i++)
        {
            GetLocalGamma(i, gamma);
            GetLocalF(LocB, i);
            MassMxMultVec(LocB);
            for (int k = 0; k < 9; k++)
            {
                _GlobB[_ElBonds[0][i][k]] += LocB[k];
            }
            GetLocG(i, LocG);
            for (int k = 0; k < 9; k++)
            {
                _di[_ElBonds[0][i][k]] += gamma * _LocMass[k][k] + LocG[k][k];
            }

            int Index;
            for (int k = 1; k < 9; k++)
            {
                for (int j = 0; j < k; j++)
                {
                    GetIndex(_ElBonds[0][i][k], _ElBonds[0][i][j], Index);
                    _al[Index] += gamma * _LocMass[k][j] + LocG[k][j];
                    _au[Index] += gamma * _LocMass[k][j] + LocG[k][j];
                }
            }
        }
    }

    //Вовзращает положение элемента в векторах al/au для его вставки в глобальную матрицу
    void GetIndex(int i, int j, int &Index)
    {
        Index = _ig[i];
        while (_jg[Index] != j)
        {
            Index++;
        }
    }

    //Собирает из локальных матриц жёсткости, полученных в ходе разложение коэффициента диффузии по билинейным базисным функциям
    void GetLocG(int num, std::vector<std::vector<double>> &LocG)
    {
        std::vector<double> LocPhiLam = std::vector<double>(4);
        GetLocLam(num, LocPhiLam);
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
    void GetLocLam(int num, std::vector<double> &LocPhiLam)
    {
        for (int i = 0; i < 4; i++)
        {
            LocPhiLam[i] = Lambda(_Grid[_ElBonds[1][num][i]][0], _Grid[_ElBonds[1][num][i]][1]);
        }
        int k = 0;
    }

    //Возвращает усреднённый коэффициент при матрице масс для элемента
    void GetLocalGamma(int num, double &gamma)
    {
        double hx, hy;
        double g = 0;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        for (int i = 0; i < 9; i++)
        {
            g += Gamma(_Grid[_ElBonds[1][num][0]][0] + hx * (i % 3), _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3));
        }
        gamma = g / 9.;
    }

    //Умножает вектор на матрицу масс для вставки в правую часть СЛАУ
    void MassMxMultVec(std::vector<double> &f)
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
    void GetLocalF(std::vector<double> &LocB, int num)
    {
        double hx, hy;
        hx = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]) / 2.;
        std::vector<double> res = std::vector<double>(9);
        for (int i = 0; i < 9; i++)
        {
            double x = _Grid[_ElBonds[1][num][0]][0] + hx * (i % 3);
            double y = _Grid[_ElBonds[1][num][0]][1] + hy * (i / 3);
            res[i] = F(x, y);
        }
        LocB = res;
    }

    //Учёт второго краевого условия
    void SecondBoundary()
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
                    GetLocalTheta(num, h / 2.0, Border, LocTheta);
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
                    GetLocalTheta(num, h / 2., Border, LocTheta);
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
                    GetLocalTheta(num, h / 2., Border, LocTheta);
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
                    GetLocalTheta(num, h / 2., Border, LocTheta);
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
    void GetLocalTheta(int num, double h, int Border, std::vector<double> &LocTheta)
    {

        for (int i = 0; i < 3; i++)
        {
            switch (Border)
            {
            case 0:
            {
                LocTheta[i] = Theta(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * i, Border);
                break;
            }
            case 1:
            {
                LocTheta[i] = Theta(_Grid[_ElBonds[1][num][0]][0] + h * i, _Grid[_ElBonds[1][num][0]][1], Border);
                break;
            }
            case 2:
            {
                LocTheta[i] = Theta(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * i, Border);
                break;
            }
            case 3:
            {
                LocTheta[i] = Theta(_Grid[_ElBonds[1][num][2]][0] + h * i, _Grid[_ElBonds[1][num][2]][1], Border);
                break;
            }
            default:
                break;
            }
        }
    }

    //Возвращает локальный вектор значений функции Убета на ребре для учёта в третем краевом условии
    void GetLocalUbeta(int num, double h, int Border, std::vector<double> &LocUbeta)
    {
        for (int i = 0; i < 3; i++)
        {
            switch (Border)
            {
            case 0:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * i, Border);
                break;
            }
            case 1:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][0]][0] + h * i, _Grid[_ElBonds[1][num][0]][1], Border);
                break;
            }
            case 2:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * i, Border);
                break;
            }
            case 3:
            {
                LocUbeta[i] = Ubeta(_Grid[_ElBonds[1][num][2]][0] + h * i, _Grid[_ElBonds[1][num][2]][1], Border);
                break;
            }
            default:
                break;
            }
        }
    }

    //Учёт третьего краевого условия
    void ThirdBoundary()
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
                    double beta = Beta();
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    InsertLocBondA(num, h, Border, LocBondA);
                    GetLocalUbeta(num, h / 2., Border, LocUbeta);
                    _GlobB[_ElBonds[0][num][0]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][3]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][6]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 1:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i;
                    double beta = Beta();
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    InsertLocBondA(num, h, Border, LocBondA);
                    GetLocalUbeta(num, h / 2., Border, LocUbeta);
                    _GlobB[_ElBonds[0][num][0]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][1]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][2]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 2:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i * (_nx - 1) + _nx - 2;
                    double beta = Beta();
                    double h = (_Grid[_ElBonds[1][num][2]][1] - _Grid[_ElBonds[1][num][0]][1]);
                    InsertLocBondA(num, h, Border, LocBondA);
                    GetLocalUbeta(num, h / 2., Border, LocUbeta);
                    _GlobB[_ElBonds[0][num][2]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][5]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            case 3:
            {
                int num;
                for (int i = 0; i < (_nx - 1); i++)
                {
                    num = i + (_nx - 1) * (_ny - 2);
                    double beta = Beta();
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]);
                    InsertLocBondA(num, h, Border, LocBondA);
                    GetLocalUbeta(num, h / 2., Border, LocUbeta);
                    _GlobB[_ElBonds[0][num][6]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][7]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
                    _GlobB[_ElBonds[0][num][8]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
                }
                break;
            }
            }
        }
    }

    //Вставляет в глобальную матрицу локальную матрицу ребра при учёте третьего краевого условия
    void InsertLocBondA(int num, double h, int Border, std::vector<std::vector<double>> LocBondA)
    {
        switch (Border)
        {
        case 0:
        {
            double beta = Beta();
            _di[_ElBonds[0][num][0]] += h * beta * LocBondA[0][0];
            _di[_ElBonds[0][num][3]] += h * beta * LocBondA[1][1];
            _di[_ElBonds[0][num][6]] += h * beta * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    GetIndex(_ElBonds[0][num][i * 3], _ElBonds[0][num][3 * j], Index);
                    _al[Index] += h * beta * LocBondA[i][j];
                    _au[Index] += h * beta * LocBondA[i][j];
                }
            }
            break;
        }
        case 1:
        {
            double beta = Beta();
            _di[_ElBonds[0][num][0]] += h * beta * LocBondA[0][0];
            _di[_ElBonds[0][num][1]] += h * beta * LocBondA[1][1];
            _di[_ElBonds[0][num][2]] += h * beta * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    GetIndex(_ElBonds[0][num][i], _ElBonds[0][num][j], Index);
                    _al[Index] += h * beta * LocBondA[i][j];
                    _au[Index] += h * beta * LocBondA[i][j];
                }
            }
            break;
        }
        case 2:
        {
            double beta = Beta();
            _di[_ElBonds[0][num][2]] += h * beta * LocBondA[0][0];
            _di[_ElBonds[0][num][5]] += h * beta * LocBondA[1][1];
            _di[_ElBonds[0][num][8]] += h * beta * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    GetIndex(_ElBonds[0][num][3 * i + 2], _ElBonds[0][num][3 * j + 2], Index);
                    _al[Index] += h * beta * LocBondA[i][j];
                    _au[Index] += h * beta * LocBondA[i][j];
                }
            }
            break;
        }
        case 3:
        {
            double beta = Beta();
            _di[_ElBonds[0][num][6]] += h * beta * LocBondA[0][0];
            _di[_ElBonds[0][num][7]] += h * beta * LocBondA[1][1];
            _di[_ElBonds[0][num][8]] += h * beta * LocBondA[2][2];
            int Index;
            for (int i = 1; i < 3; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    GetIndex(_ElBonds[0][num][6 + i], _ElBonds[0][num][6 + j], Index);
                    _al[Index] += h * beta * LocBondA[i][j];
                    _au[Index] += h * beta * LocBondA[i][j];
                }
            }
            break;
        }
        }
    }

    //Учёт первого краевого условия (Зануление строки)
    void FirstBoundary()
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
                        Ugi = Ug(_Grid[_ElBonds[1][num][0]][0], _Grid[_ElBonds[1][num][0]][1] + h * j, Border);
                        NullifyStr(_ElBonds[0][num][3 * j], Ugi);
                        _di[_ElBonds[0][num][3 * j]] = 1.;
                        _GlobB[_ElBonds[0][num][3 * j]] = _di[_ElBonds[0][num][3 * j]] * Ugi;
                    }
                }
                break;
            }
            case 1:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i;
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][0]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][0]][0] + h * j, _Grid[_ElBonds[1][num][0]][1], Border);
                        NullifyStr(_ElBonds[0][num][j], Ugi);
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
                        Ugi = Ug(_Grid[_ElBonds[1][num][1]][0], _Grid[_ElBonds[1][num][1]][1] + h * j, Border);
                        NullifyStr(_ElBonds[0][num][3 * j + 2], Ugi);
                        _di[_ElBonds[0][num][3 * j + 2]] = 1.;
                        _GlobB[_ElBonds[0][num][3 * j + 2]] = _di[_ElBonds[0][num][3 * j + 2]] * Ugi;
                    }
                }
                break;
            }
            case 3:
            {
                int num;
                for (int i = 0; i < (_ny - 1); i++)
                {
                    num = i + (_nx - 2) * (_ny - 1);
                    double h = (_Grid[_ElBonds[1][num][1]][0] - _Grid[_ElBonds[1][num][0]][1]) / 2.0;
                    double Ugi;
                    double C = 1e+30;
                    for (int j = 0; j < 3; j++)
                    {
                        Ugi = Ug(_Grid[_ElBonds[1][num][2]][0] + h * j, _Grid[_ElBonds[1][num][2]][1], Border);
                        NullifyStr(_ElBonds[0][num][6 + j], Ugi);
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
    void NullifyStr(int Node, double Ugi)
    {
        int CurNode;
        for (int i = _ig[Node]; i < _ig[Node + 1]; i++)
        {
            //_GlobB[_jg[i]] -= Ugi * _al[i];
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
                    //_GlobB[i] -= _al[_ig[i] + k] * Ugi;
                    _au[_ig[i] + k] = 0.;
                }
                k++;
                if (k < (_ig[i + 1] - _ig[i]))
                    CurNode = _jg[_ig[i] + k];
            }
        }
    }

    //Инициализирует начальное приближение (Нули)
    void InitApprox()
    {
        _Q.resize((2 * _nx - 1) * (2 * _ny - 1));
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
    std::vector<double> Mult(std::vector<double> &v)
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
    void LoS_precond()
    {
        std::vector<double> _dif = _di;
        std::vector<double> _auf = _au;
        std::vector<double> _alf = _al;
        LUsq(_dif, _alf, _auf);
        int k = 0;
        std::vector<double> buf = Mult(_Q);
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            buf[i] = _GlobB[i] - buf[i];
        }
        std::vector<double> r = LUDirect(buf, _dif, _alf);
        double error = scalar_prod(r, r);
        double error1 = error + 1;
        std::vector<double> z = LUReverse(r, _dif, _auf);
        buf = Mult(z);
        std::vector<double> p = LUDirect(buf, _dif, _alf);
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
                _Q[i] += alpha * z[i];
                r[i] -= alpha * p[i];
            }
            std::vector<double> Ur = LUReverse(r, _dif, _auf);
            buf = Mult(Ur);
            buf = LUDirect(buf, _dif, _alf);
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
    std::vector<double> LUDirect(const std::vector<double> &b, std::vector<double> &_dif, std::vector<double> &_alf)
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
    std::vector<double> LUReverse(const std::vector<double> &b, std::vector<double> &_dif, std::vector<double> &ggu_f)
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
    void print_Mx(std::vector<double> &di, std::vector<double> &ggu, std::vector<double> &ggl)
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

    //Выводит решение СЛАУ, значение искомой функции в узлах и их разницу
    void PrintError()
    {
        double hx, hy;
        hx = (_Grid[_ElBonds[1][0][1]][0] - _Grid[_ElBonds[1][0][0]][0]) / 2.;
        hy = (_Grid[_ElBonds[1][0][2]][1] - _Grid[_ElBonds[1][0][0]][1]) / 2.;
        std::cout << std::setw(17) << "(x;y)"
                  << " | " << std::setw(12) << "q(i)"
                  << " | " << std::setw(12)
                  << "u(i)"
                  << " | " << std::setw(12)
                  << "q(i) - u(i)"
                  << " | " << std::endl;
        std::cout << "__________________|______________|______________|______________|" << std::endl;
        for (int i = 0; i < (2 * _nx - 1) * (2 * _ny - 1); i++)
        {
            double x = _Grid[_ElBonds[1][0][0]][0] + hx * (i % (2 * _nx - 1));
            double y = _Grid[_ElBonds[1][0][0]][1] + hy * (i / (2 * _ny - 1));
            std::cout << "(" << std::setw(7) << x << ";"
                      << std::setw(7) << y << ")"
                                              " | "
                      << std::setw(12) << _Q[i] << " | " << std::setw(12)
                      << U(x, y) << " | " << std::setw(12) << _Q[i] - U(x, y) << " | " << std::endl;
        }
    }
};

int main()
{
    matrix MX;
    MX.makingGrid();     //Загружаем область решения и разбиваем её на конечные элементы
    MX.SetProfile();     //Инициализируем профиль глобальной матрицы
    MX.LocalBuild();     //Собираем локальную матрицу
    MX.GlobalBuild();    //Собираем глобальную матрицу
    MX.SecondBoundary(); //Учитываем второе краевое условие
    MX.ThirdBoundary();  //Учитываем третье краевое условие
    MX.FirstBoundary();  //Учитываем первое краевое условие
    MX.InitApprox();     //Инициализируем начальное приближение
    MX.LoS_precond();    //Решаем СЛАУ
    MX.PrintError();     //Выводим ответ

    return 0;
}
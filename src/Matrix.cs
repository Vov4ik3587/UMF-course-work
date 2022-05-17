using System.Drawing;
using static System.Math;
    
namespace kursach;

/// <summary>
/// Матрица в разреженно-строчном формате (CSR)
/// 
/// </summary>
public class Matrix
{
    /// <summary>
    /// Массив ненулевых элементов нижнего треугольника матрицы по строкам
    /// </summary>
    /// <returns></returns>
    public double[] Ggl;

    /// <summary>
    /// Массив ненулевых элементов верхнего треугольника матрицы по столбцам
    /// </summary>
    public double[] Ggu;

    /// Массив диагональных элементов матрицы 
    public readonly double[] Di;

    /// <summary>
    /// Массив номеров элементов матриц Ggl и Ggu, с которых начинается i-ая строка нижнего треугольника
    /// и i-ый стобец верхнего треугольника матрицы
    /// Длина массива N+1, где N - размерность матрицы А
    /// </summary>
    public int[] Ig;

    /// <summary>
    /// 
    /// </summary>
    public int[] Jg;

    public Matrix()
    {
        Ggl = default!;
        Ggu = default!;
        Di = default!;
        Ig = default!;
        Jg = default!;
        Decomposed = default!;
        Size = default!;
    }

    public Matrix(double[] ggl, double[] ggu, double[] di, int[] ig, int[] jg, int size, bool decomposed)
    {
        Ggl = ggl;
        Ggu = ggu;
        Di = di;
        Ig = ig;
        Jg = jg;
        Size = size;
        Decomposed = decomposed;
    }

    /// <summary>
    /// Было LU-разложение или нет  
    /// </summary>
    public bool Decomposed { get; private set; }

    /// <summary>
    /// Размерность матрицы
    /// </summary>
    public int Size { get; }

    /// <summary>
    /// LU-разложение матрицы. Она изменяет сам объект, а не создает его копию
    /// </summary>
    /// <exception cref="DivideByZeroException"> If diagonal element is zero </exception>
    public void Factorize()
    {
        for (var i = 0; i < Size; i++)
        {
            var sumDi = 0.0;

            var i0 = Ig[i];
            var i1 = Ig[i + 1];

            for (var k = i0; k < i1; k++)
            {
                var j = Jg[k];
                var j0 = Ig[j];
                var j1 = Ig[j + 1];

                var iK = i0;
                var kJ = j0;

                var sumL = 0.0;
                var sumU = 0.0;

                while (iK < k && kJ < j1)
                {
                    if (Jg[iK] == Jg[kJ])
                    {
                        sumL += Ggl[iK] * Ggu[kJ];
                        sumU += Ggu[iK] * Ggl[kJ];
                        iK++;
                        kJ++;
                    }
                    else
                    {
                        if (Jg[iK] > Jg[kJ])
                        {
                            kJ++;
                        }
                        else
                        {
                            iK++;
                        }
                    }
                }

                if (Di[j] == 0.0)
                {
                    throw new DivideByZeroException($"Di[{j}] has thrown at pos {i} {j}");
                }
                Ggl[k] = (Ggl[k] - sumL) / Di[j];
                Ggu[k] = (Ggu[k] - sumU) / Di[j];
                
                sumDi += Ggl[k] * Ggu[k];
            }

            Di[i] = Sqrt(Di[i] - sumDi);
        }

        Decomposed = true;
    }
}


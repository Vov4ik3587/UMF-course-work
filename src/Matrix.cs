namespace kursach;

public class Matrix
{
    public double[] Aelem;  // Все ненулевые элементы матрицы по строкам
    public double[] Jptr; // Номера столбцов соответствующих ненулевых элементов
    public double[] Iptr; // Всем понятный массив ia, не знаю, как его назвать

    public int Size { get; }
    public Matrix()
    {
        Aelem = default!;
        Jptr = default!;
        Iptr = default!;
    }

    public Matrix(double[] aelem, double[] jptr, double[] iptr, int size)
    {
        Aelem = aelem ?? throw new ArgumentNullException(nameof(aelem));
        Jptr = jptr ?? throw new ArgumentNullException(nameof(jptr));
        Iptr = iptr ?? throw new ArgumentNullException(nameof(iptr));
        Size = size;
    }

}
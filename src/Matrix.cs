namespace kursach;


/// <summary>
/// Матрица в разреженно-строчном формате
/// </summary>
public class Matrix
{
    public double[] Aelem { get; }
    public double[] Jptr { get; } // Номера столбцов соответствующих ненулевых элементов
    public double[] Ia { get; } 
    
}

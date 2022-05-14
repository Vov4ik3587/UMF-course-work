namespace kursach;

public class BasisFunc
{
    public static Func<double, double, double>[] Psi =
    {
        (x,y) => x * y,

    };

    public static Func<double, double, double>[] Lambda =
    {
        (x, y) => x * y,
        (x, y) => (1 - x) * y,
        (x, y) => x * (1 - y),
        (x, y) => (1 - x) * (1 - y),
    };
}
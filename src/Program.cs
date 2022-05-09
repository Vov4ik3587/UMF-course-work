
namespace kursach
{
    internal static class Program
    {
        /// <summary>
        /// Запускает выбранный вариант курсача: станционарную или нестанционарную краевую задачу
        /// </summary>
        /// <param name="args">
        /// 5-semester - станционарная задача
        /// 6-semester - нестанционарная задача
        /// </param>
        public static void Main(string[] args)
        {
            switch (args[0])
            {
                case "5-semester":
                    Console.WriteLine("Решение станционарной задачи из 5 семестра начинается...");
                    
                    break;
                case "6-semester": 
                    Console.WriteLine("Решение нестанционарной задачи из 6 семестра начинается...");
                    break;
            }
        }
    }
}
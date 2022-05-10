
using System.Text.Json;
using System.Text.Json.Nodes;

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
                    StationaryTask();
                    break;
                case "6-semester": 
                    Console.WriteLine("Решение нестанционарной задачи из 6 семестра начинается...");
                    NoStationaryTask();
                    break;
            }
        }

        public static void StationaryTask()
        {
            Console.WriteLine("Начинаем считывание данных...");
            var area = JsonSerializer.Deserialize<Area>(File.ReadAllText("Input-data/area.json"));
        }
        public static void NoStationaryTask()
        {
            Console.WriteLine("Реализации еще нет...Будет позже");
        }
    }
}
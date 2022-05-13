namespace kursach
{
    public class All_Math
    {
        /// <summary>
        /// Matrix multiplication by vector
        /// </summary>
        /// <param name="m"> Matrix, stored in 2-dim double array</param>
        /// <param name="b"> Vector, that multiplies matrix </param>
        /// <returns> Resulting vector of multiplication </returns>
        public static double[] MatMul(double[,] m, double[] b)
        {
            var columns = m.GetLength(0);
            var rows = m.GetLength(1);

            var res = new double[rows];

            for (var i = 0; i < rows; i++)
            {
                for (var j = 0; j < columns; j++)
                {
                    res[i] += m[i, j] * b[j];
                }
            }

            return res;
        }

        /// <summary>
        /// Matrix multiplication by vector
        /// </summary>
        /// <param name="m"> Matrix, stored in sparse format</param>
        /// <param name="b"> Vector, that multiplies matrix </param>
        /// <returns> Resulting vector of multiplication </returns>
        public static double[] MatMul(Matrix m, double[] b)
        {
            if (m.Size != b.Length)
            {
                throw new Exception($"[ERR] Different sizes. Matrix size = {m.Size}, vector size = {b.Length}");
            }

            var res = new double[b.Length];

            for (var i = 0; i < b.Length; i++)
            {
                res[i] = m.Di[i] * b[i];

                for (var j = m.Ig[i]; j < m.Ig[i + 1]; j++)
                {
                    res[i] += m.Ggl[j] * b[m.Jg[j]];
                    res[m.Jg[j]] += m.Ggu[j] * b[i];
                }
            }

            return res;
        }
        // Добавить Сложение, произведение матриц в разреженном формате

        // TODO: Страница кирпича 217 
                    // Вычисление мю(i)
                    // Вычисление ню(i)
    }
}
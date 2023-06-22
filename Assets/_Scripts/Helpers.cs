using System.Collections.Generic;
using System.Linq;

using MathNet.Numerics.LinearAlgebra;
public static class Helpers
{
    public static int MultiplyArray(IEnumerable<int> numbers)
    {
        return numbers.Aggregate(1, (current, value) => current * value);
    }
    
    public static double GetDeterminant(List<List<double>> matrix)
    {
        double determinant = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                             matrix[0][1] * matrix[1][2] * matrix[2][0] +
                             matrix[1][0] * matrix[0][2] * matrix[2][1] -
                             matrix[2][0] * matrix[1][1] * matrix[0][2] -
                             matrix[0][1] * matrix[1][0] * matrix[2][2] -
                             matrix[0][0] * matrix[1][2] * matrix[2][1];
    
        return determinant;
    }
    
    public static Matrix<double> CreateMatrixFrom2DList(List<List<double>> list)
    {
        int rows = list.Count;
        int columns = list[0].Count;
        double[,] array = new double[rows, columns];

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                array[i, j] = list[i][j];
            }
        }

        return Matrix<double>.Build.DenseOfArray(array);
    }

    public static MathNet.Numerics.LinearAlgebra.Vector<double> CreateVectorFromList(List<double> list)
    {
        return MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(list.ToArray());
    }

    public static List<double> CreateListFromVector(MathNet.Numerics.LinearAlgebra.Vector<double> vector)
    {
        return vector.ToList();
    }
}
using System;
using System.Collections.Generic;

public static class GlobalMatrix
{
    public static void InitGlobals(List<List<double>> globalMatrix, List<double> vecF, List<List<List<double>>> allMGE, List<List<double>> allFE, List<List<int>> ZP, List<List<int>> NT)
    {
        for (int q = 0; q < allMGE.Count; q++)
        {
            for (int i = 0; i < allMGE[0].Count; i++)
            {
                for (int j = 0; j < allMGE[0].Count; j++)
                {
                    int compx = i / 20;
                    int compy = j / 20;
                    int i1 = NT[i % 20][q] * 3 + compx;
                    int i2 = NT[j % 20][q] * 3 + compy;
                    globalMatrix[i1][i2] += allMGE[q][i][j];
                }
            }
        }

        for (int q = 0; q < ZP.Count; q++)
        {
            for (int i = 0; i < allFE[0].Count; i++)
            {
                int comp = i / 20;
                int i1 = NT[i % 20][ZP[q][0]] * 3 + comp;
                vecF[i1] += allFE[q][i];
            }
        }
    }

    public static void FixateSides(List<List<double>> globalMatrix, List<List<int>> ZU, List<List<int>> NT)
    {
        for (int i = 0; i < ZU.Count; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                int index = NT[Constants.SideToPoint[ZU[i][1]][j]][ZU[i][0]] * 3;

                globalMatrix[index][index] = Math.Pow(10, 30);
                globalMatrix[index + 1][index + 1] = Math.Pow(10, 30);
                globalMatrix[index + 2][index + 2] = Math.Pow(10, 30);
            }
        }
    }
}
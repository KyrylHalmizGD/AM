using System;
using System.Collections.Generic;
using System.Linq;


public static class ForceVector
{
    public static List<List<List<double>>> CreateDPSIET()
    {
        List<List<double>> nineNodes = GaussNodes.GetNineGaussianNodes();
        List<List<List<double>>> DPSIET = new List<List<List<double>>>();

        for (int i = 0; i < 9; i++)
        {
            List<List<double>> row = new List<List<double>>();
            for (int j = 0; j < 2; j++)
            {
                List<double> col = new List<double>();
                for (int k = 0; k < 8; k++)
                {
                    double value = 0.0;
                    if (j == 0)
                    {
                        if (k < 4)
                            value = DeltaPsi.DeltaPsiEta4(nineNodes[i][0], nineNodes[i][1], k);
                        else if (k == 4)
                            value = DeltaPsi.DeltaPsiEta57(nineNodes[i][0], nineNodes[i][1], 4);
                        else if (k == 6)
                            value = DeltaPsi.DeltaPsiEta57(nineNodes[i][0], nineNodes[i][1], 6);
                        else if (k == 5)
                            value = DeltaPsi.DeltaPsiEta68(nineNodes[i][0], nineNodes[i][1], 5);
                        else if (k == 7)
                            value = DeltaPsi.DeltaPsiEta68(nineNodes[i][0], nineNodes[i][1], 7);
                    }
                    else if (j == 1)
                    {
                        if (k < 4)
                            value = DeltaPsi.DeltaPsiTau4(nineNodes[i][0], nineNodes[i][1], k);
                        else if (k == 4)
                            value = DeltaPsi.DeltaPsiTau57(nineNodes[i][0], nineNodes[i][1], 4);
                        else if (k == 6)
                            value = DeltaPsi.DeltaPsiTau57(nineNodes[i][0], nineNodes[i][1], 6);
                        else if (k == 5)
                            value = DeltaPsi.DeltaPsiTau68(nineNodes[i][0], nineNodes[i][1], 5);
                        else if (k == 7)
                            value = DeltaPsi.DeltaPsiTau68(nineNodes[i][0], nineNodes[i][1], 7);
                    }
                    col.Add(value);
                }
                row.Add(col);
            }
            DPSIET.Add(row);
        }

        return DPSIET;
    }   
    
    public static List<List<List<double>>> CreateDPSIXYZ(List<List<List<double>>> DPSIET, List<List<double>> AKT, List<List<int>> NT, int finiteElementIndex, int side)
    {
        List<List<List<double>>> DPSIXYZ = new List<List<List<double>>>();

        for (int i = 0; i < 9; i++)
        {
            List<List<double>> row = new List<List<double>>();
            for (int j = 0; j < 2; j++)
            {
                List<double> col = new List<double>();
                for (int k = 0; k < 3; k++)
                {
                    double sum = 0.0;
                    for (int s = 0; s < 8; s++)
                    {
                        int pointIndex = Constants.SideToPoint[side][s];
                        int ntValue = NT[pointIndex][finiteElementIndex];
                        double aktValue = AKT[k][ntValue];
                        double dpsietValue = DPSIET[i][j][s];
                        sum += aktValue * dpsietValue;
                    }
                    col.Add(sum);
                }
                row.Add(col);
            }
            DPSIXYZ.Add(row);
        }

        return DPSIXYZ;
    }

public static List<double> CreateFe1(List<List<List<double>>> DPSIET, List<List<double>> AKT, List<List<int>> NT, List<int> ZPi, Dictionary<Tuple<double, double>, int> permutationIndex)
{
    double Pn = ZPi[2];
    List<double> fe1 = Enumerable.Repeat(0.0, 8).ToList();
    List<List<List<double>>> DPSIXYZ = CreateDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1]);
    for (int i = 0; i < 8; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            for (int n = 0; n < 3; n++)
            {
                Tuple<double, double> key = Tuple.Create(Constants.Xg[m], Constants.Xg[n]);
                int pIndex = permutationIndex[key];
                double nul = DPSIXYZ[pIndex][0][1] * DPSIXYZ[pIndex][1][2] - DPSIXYZ[pIndex][0][2] * DPSIXYZ[pIndex][1][1];
                fe1[i] += Constants.C[m] * Constants.C[n] * Pn * nul * Psi.Psii(Constants.Xg[m], Constants.Xg[n], i, i);
            }
        }
    }
    return fe1;
}

public static List<double> CreateFe2(List<List<List<double>>> DPSIET, List<List<double>> AKT, List<List<int>> NT, List<int> ZPi, Dictionary<Tuple<double, double>, int> permutationIndex)
{
    double Pn = ZPi[2];
    List<double> fe2 = Enumerable.Repeat(0.0, 8).ToList();
    List<List<List<double>>> DPSIXYZ = CreateDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1]);
    for (int i = 0; i < 8; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            for (int n = 0; n < 3; n++)
            {
                Tuple<double, double> key = Tuple.Create(Constants.Xg[m], Constants.Xg[n]);
                int pIndex = permutationIndex[key];
                fe2[i] += Constants.C[m] * Constants.C[n] * Pn * (DPSIXYZ[pIndex][0][2] * DPSIXYZ[pIndex][1][0] - DPSIXYZ[pIndex][0][0] * DPSIXYZ[pIndex][1][2]) * Psi.Psii(Constants.Xg[m], Constants.Xg[n], i, i);
            }
        }
    }
    return fe2;
}

public static List<double> CreateFe3(List<List<List<double>>> DPSIET, List<List<double>> AKT, List<List<int>> NT, List<int> ZPi, Dictionary<Tuple<double, double>, int> permutationIndex)
{
    double Pn = ZPi[2];
    List<double> fe3 = Enumerable.Repeat(0.0, 8).ToList();
    List<List<List<double>>> DPSIXYZ = CreateDPSIXYZ(DPSIET, AKT, NT, ZPi[0], ZPi[1]);
    for (int i = 0; i < 8; i++)
    {
        for (int m = 0; m < 3; m++)
        {
            for (int n = 0; n < 3; n++)
            {
                Tuple<double, double> key = Tuple.Create(Constants.Xg[m], Constants.Xg[n]);
                int pIndex = permutationIndex[key];
                fe3[i] += Constants.C[m] * Constants.C[n] * Pn * (DPSIXYZ[pIndex][0][0] * DPSIXYZ[pIndex][1][1] - DPSIXYZ[pIndex][0][1] * DPSIXYZ[pIndex][1][0]) * Psi.Psii(Constants.Xg[m], Constants.Xg[n], i, i);
            }
        }
    }
    return fe3;
}

public static List<double> CreateFE(List<List<List<double>>> DPSIET, List<List<double>> AKT, List<List<int>> NT, List<int> ZPi, Dictionary<Tuple<double, double>, int> permutationIndex)
{
    List<double> fe1 = CreateFe1(DPSIET, AKT, NT, ZPi, permutationIndex);
    List<double> fe2 = CreateFe2(DPSIET, AKT, NT, ZPi, permutationIndex);
    List<double> fe3 = CreateFe3(DPSIET, AKT, NT, ZPi, permutationIndex);
    List<double> fe = Enumerable.Repeat(0.0, 60).ToList();
    for (int i = 0; i < 8; i++)
    {
        fe[Constants.SideToPoint[ZPi[1]][i]] = fe1[i];
        fe[Constants.SideToPoint[ZPi[1]][i] + 20] = fe2[i];
        fe[Constants.SideToPoint[ZPi[1]][i] + 40] = fe3[i];
    }
    return fe;
}

}
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


internal abstract class Program
{
    private static void Main()
    {
        var xyzLen = new List<int>() { 4, 2, 4 };// довжина початкового цилідру 
        var xyzQuan = new List<int>() { 2, 1, 2 };// кількість розбиттів по x,y,z відповідно 
        
        List<List<int>> zp = new List<List<int>> { new() { 2, 6, -1 }, new() { 3, 6, -1 } };// визначення сили 

        List<List<int>> zu = new List<List<int>> { new() { 0, 5 }, new() { 1, 5 } };// закріплені елементи
        var npq = Helpers.MultiplyArray(xyzQuan);
        
        Console.WriteLine("npq = " + npq);

        List<List<double>> akt = Net.GenerateAkt(xyzLen,xyzQuan);

        List<List<int>> nt = Net.GenerateNt(akt, xyzQuan);
        

        List<List<List<double>>> dfiAbg = Stiffness.CreateDfiAbg();

        List<List<List<List<double>>>> jacobiansForElement = Stiffness.GenerateJacobiansForElement(nt,akt,dfiAbg,npq);

        List<List<double>> jacobianValuesForElement = Stiffness.GetJacobianValuesForElement(jacobiansForElement,npq);

        List<List<List<List<double>>>> dfiXyz = Stiffness.CreateDfiXyz(jacobiansForElement, dfiAbg, npq);
        
        List<List<double>> permutations = GaussNodes.GetTwentySevenGaussianNodes();
        
        Dictionary<(double, double, double), int> permutationIndex = new Dictionary<(double, double, double), int>();
        for (int i = 0; i < 27; i++)
        {
            permutationIndex.Add((permutations[i][0], permutations[i][1], permutations[i][2]), i);
        }

        List<List<List<double>>> allAe11 = Stiffness.CreateAe11(dfiXyz, jacobianValuesForElement, npq, permutationIndex);
        List<List<List<double>>> allAe22 = Stiffness.CreateAe22(dfiXyz, jacobianValuesForElement, npq, permutationIndex);
        List<List<List<double>>> allAe33 = Stiffness.CreateAe33(dfiXyz, jacobianValuesForElement, npq, permutationIndex);
        List<List<List<double>>> allAe12 = Stiffness.CreateAe12(dfiXyz, jacobianValuesForElement, npq, permutationIndex);
        List<List<List<double>>> allAe13 = Stiffness.CreateAe13(dfiXyz, jacobianValuesForElement, npq, permutationIndex);
        List<List<List<double>>> allAe23 = Stiffness.CreateAe23(dfiXyz, jacobianValuesForElement, npq, permutationIndex);

        var allMge = Stiffness.CreateMGE(allAe11, allAe22, allAe33, allAe12, allAe13, allAe23, npq);

        var dpsiet = ForceVector.CreateDPSIET();
        
        List<List<double>> permutations2d = GaussNodes.GetNineGaussianNodes();
        Dictionary<Tuple<double, double>, int> permutationIndex2d = new Dictionary<Tuple<double, double>, int>();

        for (int i = 0; i < 9; i++)
        {
            permutationIndex2d[Tuple.Create(permutations2d[i][0], permutations2d[i][1])] = i;
        }

        List<List<double>> allFe = zp.Select(t => ForceVector.CreateFE(dpsiet, akt, nt, t, permutationIndex2d)).ToList();
        
        // Making all MGE full matrices (tmp to be usable with current gaussian elimination)
        foreach (var t in allMge)
        {
            for (int i = 0; i < allMge[0].Count; i++)
            {
                for (int j = i; j < allMge[0][0].Count; j++)
                {
                    t[j][i] = t[i][j];
                }
            }
        }
        
        List<List<double>> globalMG = new List<List<double>>();
        var nqp = (akt[0].Count);
        List<double> globalF = new List<double>();

       
        for (int i = 0; i < nqp * 3; i++)
        {
            List<double> row = new List<double>();
            for (int j = 0; j < nqp * 3; j++)
            {
                row.Add(0.0);
            }
            globalMG.Add(row);
        }
        
        // Initialize globalF with zeros
        for (int i = 0; i < nqp * 3; i++)
        {
            globalF.Add(0);
        }

        GlobalMatrix.InitGlobals(globalMG, globalF, allMge, allFe, zp, nt);
        GlobalMatrix.FixateSides(globalMG, zu, nt);

        
        Vector<double> finalRes = DenseVector.Create(globalF.Count, 0);

        Matrix<double> globalMGMatrix = DenseMatrix.OfRows(globalMG);
        Vector<double> globalFVector = DenseVector.OfEnumerable(globalF);

        finalRes = globalMGMatrix.Solve(globalFVector);


        var coef = 30;
        
        for (var i = 0; i < akt[0].Count; i++)
        {
            akt[0][i] += coef * finalRes[i * 3];
            akt[1][i] += coef * finalRes[i * 3 + 1];
            akt[2][i] += coef * finalRes[i * 3 + 2];
        }
        
        Console.WriteLine("success");
    }
}
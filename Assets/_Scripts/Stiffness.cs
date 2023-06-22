using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
public static class Stiffness
{
    public static List<List<List<double>>> CreateDfiAbg() // 27*3*20 27 вузлів гауса * 3(альфа, бета, гамма) * 20 точок(функцій)
    {
        var permutations = GaussNodes.GetTwentySevenGaussianNodes(); // гаусіанівські ноди 
        var dfiAbg = new List<List<List<double>>>();
        for (var i = 0; i < 27; i++)
        {
            var phiAlphaBetaGamma = new List<List<double>>();
            for (var j = 0; j < 3; j++)
            {
                var phi = new List<double>();
                for (var k = 0; k < 20; k++)
                {
                    phi.Add(0.0);
                }
                phiAlphaBetaGamma.Add(phi);
            }
            dfiAbg.Add(phiAlphaBetaGamma);
        } // заповнюємо 0ями 

        for (var i = 0; i < 27; i++)
        {
            for (var j = 0; j < 8; j++) 
            { 
                dfiAbg[i][0][j] = DeltaPhi.DeltaPhiAlpha_8(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }

            for (var j = 8; j < 20; j++)
            {
                dfiAbg[i][0][j] = DeltaPhi.DeltaPhiAlpha_12(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }

            for (var j = 0; j < 8; j++)
            {
                dfiAbg[i][1][j] = DeltaPhi.DeltaPhiBeta_8(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }

            for (var j = 8; j < 20; j++)
            {
                dfiAbg[i][1][j] = DeltaPhi.DeltaPhiBeta_12(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }

            for (var j = 0; j < 8; j++)
            {
                dfiAbg[i][2][j] = DeltaPhi.DeltaPhiGamma_8(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }

            for (var j = 8; j < 20; j++)
            {
                dfiAbg[i][2][j] = DeltaPhi.DeltaPhiGamma_12(permutations[i][0], permutations[i][1], permutations[i][2], j);
            }
        }

        return dfiAbg;
    }
    
    public static List<List<List<double>>> CreateJacobian(List<List<int>> NT, List<List<double>> AKT, List<List<List<double>>> DFIABG, int element)
    {
        int numNodes = DFIABG.Count;
        List<List<List<double>>> jacobian_matrices_27 = new List<List<List<double>>>(27);

        for (int node = 0; node < numNodes; node++)
        {
            List<List<double>> jacobian_matrix = new List<List<double>>(3);

            for (int der = 0; der < 3; der++)
            {
                List<double> jacobian_row = new List<double>(3);

                for (int coord = 0; coord < 3; coord++)
                {
                    double sum = 0.0;

                    for (int i = 0; i < 20; i++)
                    {
                        sum += AKT[coord][NT[i][element]] * DFIABG[node][der][i];
                    }

                    jacobian_row.Add(sum);
                }

                jacobian_matrix.Add(jacobian_row);
            }

            jacobian_matrices_27.Add(jacobian_matrix);
        }

        return jacobian_matrices_27;
    }
    public static List<List<List<List<double>>>> GenerateJacobiansForElement(List<List<int>> nt, List<List<double>> akt, List<List<List<double>>> dfiAbg, int npq)
    {
        List<List<List<List<double>>>> jacobiansForElement = new List<List<List<List<double>>>>(npq);
        for (int i = 0; i < npq; i++)
        {
            jacobiansForElement.Add(Stiffness.CreateJacobian(nt, akt, dfiAbg, i));
        }
        return jacobiansForElement;
    }

    public static List<List<double>> GetJacobianValuesForElement(List<List<List<List<double>>>> jacobiansForElement, int npq)
    {
        List<List<double>> jacobianValuesForElement = new List<List<double>>(npq);
        for (int i = 0; i < npq; i++)
        {
            List<double> jacobianValues = new List<double>(27);
            for (int j = 0; j < 27; j++)
            {
                double determinant = Helpers.GetDeterminant(jacobiansForElement[i][j]);
                jacobianValues.Add(determinant);
            }
            jacobianValuesForElement.Add(jacobianValues);
        }
        return jacobianValuesForElement;
    }
    
    public static List<List<List<List<double>>>> CreateDfiXyz(List<List<List<List<double>>>> jacobiansForElement, List<List<List<double>>> dfiAbg, int finiteElementsQuantity)
    {
        var dfiXyz = new List<List<List<List<double>>>>(finiteElementsQuantity);

        for (var q = 0; q < finiteElementsQuantity; q++)
        {
            var elementDfiXyz = new List<List<List<double>>>(27);

            for (int i = 0; i < 27; i++)
            {
                var nodeDfiXyz = new List<List<double>>(20);

                for (int j = 0; j < 20; j++)
                {
                    var matrix = Helpers.CreateMatrixFrom2DList(jacobiansForElement[q][i]);
                    var vector = Helpers.CreateVectorFromList(new List<double> { dfiAbg[i][0][j], dfiAbg[i][1][j], dfiAbg[i][2][j] });
                    var solution = matrix.Solve(vector);
                    nodeDfiXyz.Add(Helpers.CreateListFromVector(solution));
                }

                elementDfiXyz.Add(nodeDfiXyz);
            }

            dfiXyz.Add(elementDfiXyz);
        }

        return dfiXyz;
    }

    public static List<List<List<double>>> CreateAe11(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        var ae11 = new List<List<List<double>>>(finiteElementsQuantity);

        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            var elementAe11 = new List<List<double>>(20);

            for (int i = 0; i < 20; i++)
            {
                var nodeAe11 = new List<double>(20);

                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;

                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                (double, double, double) permutation = (Constants.Xg[m], Constants.Xg[n], Constants.Xg[k]);
                                int pIndex = permutationIndex[permutation];

                                double term = Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] +
                                    mu * (DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1] + DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2])
                                ) * jacobianValuesForElement[q][pIndex];

                                sum += term;
                            }
                        }
                    }

                    nodeAe11.Add(sum);
                }

                elementAe11.Add(nodeAe11);
            }

            ae11.Add(elementAe11);
        }

        return ae11;
    }

    public static List<List<List<double>>> CreateAe22(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        var ae22 = new List<List<List<double>>>(finiteElementsQuantity);

        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            var elementAe22 = new List<List<double>>(20);

            for (int i = 0; i < 20; i++)
            {
                var nodeAe22 = new List<double>(20);

                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;

                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                (double, double, double) permutation = (Constants.Xg[m], Constants.Xg[n], Constants.Xg[k]);
                                int pIndex = permutationIndex[permutation];

                                double term = Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1] +
                                    mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] + DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2])
                                ) * jacobianValuesForElement[q][pIndex];

                                sum += term;
                            }
                        }
                    }

                    nodeAe22.Add(sum);
                }

                elementAe22.Add(nodeAe22);
            }

            ae22.Add(elementAe22);
        }

        return ae22;
    }

    public static List<List<List<double>>> CreateAe33(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        var ae33 = new List<List<List<double>>>(finiteElementsQuantity);

        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            var elementAe33 = new List<List<double>>(20);

            for (int i = 0; i < 20; i++)
            {
                var nodeAe33 = new List<double>(20);

                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;

                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                (double, double, double) permutation = (Constants.Xg[m], Constants.Xg[n], Constants.Xg[k]);
                                int pIndex = permutationIndex[permutation];

                                double term = Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    lambda_val * (1 - nu) * DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][2] +
                                    mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][0] + DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][1])
                                ) * jacobianValuesForElement[q][pIndex];

                                sum += term;
                            }
                        }
                    }

                    nodeAe33.Add(sum);
                }

                elementAe33.Add(nodeAe33);
            }

            ae33.Add(elementAe33);
        }

        return ae33;
    }

    public static List<List<List<double>>> CreateAe12(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        var ae12 = new List<List<List<double>>>(finiteElementsQuantity);

        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            var elementAe12 = new List<List<double>>(20);

            for (int i = 0; i < 20; i++)
            {
                var nodeAe12 = new List<double>(20);

                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;

                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                (double, double, double) permutation = (Constants.Xg[m], Constants.Xg[n], Constants.Xg[k]);
                                int pIndex = permutationIndex[permutation];

                                double term = Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][1] + DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][0])
                                ) * jacobianValuesForElement[q][pIndex];

                                sum += term;
                            }
                        }
                    }

                    nodeAe12.Add(sum);
                }

                elementAe12.Add(nodeAe12);
            }

            ae12.Add(elementAe12);
        }

        return ae12;
    }

    public static List<List<List<double>>> CreateAe13(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        var ae13 = new List<List<List<double>>>(finiteElementsQuantity);

        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            var elementAe13 = new List<List<double>>(20);

            for (int i = 0; i < 20; i++)
            {
                var nodeAe13 = new List<double>(20);

                for (int j = 0; j < 20; j++)
                {
                    double sum = 0;

                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                (double, double, double) permutation = (Constants.Xg[m], Constants.Xg[n], Constants.Xg[k]);
                                int pIndex = permutationIndex[permutation];

                                double term = Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    mu * (DFIXYZ[q][pIndex][i][0] * DFIXYZ[q][pIndex][j][2] + DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][0])
                                ) * jacobianValuesForElement[q][pIndex];

                                sum += term;
                            }
                        }
                    }

                    nodeAe13.Add(sum);
                }

                elementAe13.Add(nodeAe13);
            }

            ae13.Add(elementAe13);
        }

        return ae13;
    }

    public static List<List<List<double>>> CreateAe23(List<List<List<List<double>>>> DFIXYZ, List<List<double>> jacobianValuesForElement, int finiteElementsQuantity, Dictionary<(double, double, double), int> permutationIndex)
    {
        double e = 100;
        double nu = 0.3;
        double lambda_val = e / ((1 + nu) * (1 - 2 * nu));
        double mu = e / (2 * (1 + nu));
        List<List<List<double>>> ae23 = new List<List<List<double>>>();
        for (int q = 0; q < finiteElementsQuantity; q++)
        {
            ae23.Add(new List<List<double>>());
            for (int i = 0; i < 20; i++)
            {
                ae23[q].Add(new List<double>());
                for (int j = 0; j < 20; j++)
                {
                    ae23[q][i].Add(0);
                    for (int m = 0; m < 3; m++)
                    {
                        for (int n = 0; n < 3; n++)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                int pIndex = permutationIndex[(Constants.Xg[m], Constants.Xg[n], Constants.Xg[k])];
                                ae23[q][i][j] += Constants.C[m] * Constants.C[n] * Constants.C[k] * (
                                    lambda_val * nu * DFIXYZ[q][pIndex][i][1] * DFIXYZ[q][pIndex][j][2] +
                                    mu * DFIXYZ[q][pIndex][i][2] * DFIXYZ[q][pIndex][j][1]
                                ) * jacobianValuesForElement[q][pIndex];
                            }
                        }
                    }
                }
            }
        }
        return ae23;
    }

    public static List<List<List<double>>> CreateMGE(List<List<List<double>>> allAe11, List<List<List<double>>> allAe22,
                                                     List<List<List<double>>> allAe33, List<List<List<double>>> allAe12,
                                                     List<List<List<double>>> allAe13, List<List<List<double>>> allAe23, int npq)
    {
        List<List<List<double>>> MGE = new List<List<List<double>>>();

        // Initialize MGE with zeros
        for (int q = 0; q < npq; q++)
        {
            List<List<double>> matrix = new List<List<double>>();
            for (int i = 0; i < 60; i++)
            {
                List<double> row = new List<double>();
                for (int j = 0; j < 60; j++)
                {
                    row.Add(0.0);
                }
                matrix.Add(row);
            }
            MGE.Add(matrix);
        }

        // Copy allAe11 to MGE[:, :20, :20]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    MGE[q][i][j] = allAe11[q][i][j];
                }
            }
        }

        // Copy allAe22 to MGE[:, 20:40, 20:40]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 20; i < 40; i++)
            {
                for (int j = 20; j < 40; j++)
                {
                    MGE[q][i][j] = allAe22[q][i - 20][j - 20];
                }
            }
        }

        // Copy allAe33 to MGE[:, 40:, 40:]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 40; i < 60; i++)
            {
                for (int j = 40; j < 60; j++)
                {
                    MGE[q][i][j] = allAe33[q][i - 40][j - 40];
                }
            }
        }

        // Copy allAe12 to MGE[:, :20, 20:40]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 20; j < 40; j++)
                {
                    MGE[q][i][j] = allAe12[q][i][j - 20];
                }
            }
        }

        // Copy allAe13 to MGE[:, :20, 40:]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 0; i < 20; i++)
            {
                for (int j = 40; j < 60; j++)
                {
                    MGE[q][i][j] = allAe13[q][i][j - 40];
                }
            }
        }

        // Copy allAe23 to MGE[:, 20:40, 40:]
        for (int q = 0; q < npq; q++)
        {
            for (int i = 20; i < 40; i++)
            {
                for (int j = 40; j < 60; j++)
                {
                    MGE[q][i][j] = allAe23[q][i - 20][j - 40];
                }
            }
        }

        return MGE;
    }



}
using System.Collections.Generic;
using System.Linq;

public static class Net
{
    public static List<List<double>> GenerateAkt(List<int> xyzLength, List<int> xyzQuantity) //3*npq => x,y,z координати для точок 
    {
        var xyzDistance = new List<double> {
            (double)xyzLength[0] / xyzQuantity[0],
            (double)xyzLength[1] / xyzQuantity[1],
            (double)xyzLength[2] / xyzQuantity[2]
        }; // знаходимо відстань між сусідніми точками по xyz (відстань відрізку) 

        var result = new List<List<double>> { new(), new(), new() }; // створюємо АКТ

        var layersQuantity = 1 + 2 * xyzQuantity[2]; // рахуємо кількість шарів, 
        xyzDistance[2] = xyzLength[2] / (layersQuantity - 1); // 5 точок по z - 5 шарів, 4 відстані 
        xyzDistance[1] = xyzLength[1] / ((1 + 2 * xyzQuantity[1]) - 1); // 1 шар по у  = 3 точки => 1+2 -1 = 2 відстані 
        xyzDistance[0] = xyzLength[0] / ((1 + 2 * xyzQuantity[0]) - 1); // 2 шари по х = 5 точок => 1+2*2-1 = 4 відстані

        for (var layer = 0; layer < layersQuantity; layer++)// цикл по шарах вгору
        {
            if (layer % 2 == 0) // повний шар
            {
                for (var row = 0; row < 1 + 2 * xyzQuantity[1]; row++)// проходимось по рядках 
                {
                    if (row % 2 == 0) // "повні" рядки. маємо 1 + 2 * xQuan точок, просто добавляємо величини
                    {
                        for (var i = 0; i < 1 + 2 * xyzQuantity[0]; i++)
                        {
                            result[0].Add(i * xyzDistance[0]);
                            result[1].Add(row * xyzDistance[1]);
                            result[2].Add(layer * xyzDistance[2]);
                        }
                    }
                    else // рядки на середині "половинчасті". маємо xQuan +1 точок, добавляємо величини, для х *2, оскільки немає точок на середині граней
                    {
                        for (var i = 0; i < 1 + xyzQuantity[0]; i++)
                        {
                            result[0].Add(i * xyzDistance[0] * 2);
                            result[1].Add(row * xyzDistance[1]);
                            result[2].Add(layer * xyzDistance[2]);
                        }
                    }
                }
            }
            else // "половинчатий" шар
            {
                for (var row = 0; row < 1 + xyzQuantity[1]; row++)// просто для кількості рядів добавляємо точки по х у множачи на 2, оскліьки немає точок на половинах 
                {
                    for (var i = 0; i < 1 + xyzQuantity[0]; i++)
                    {
                        result[0].Add(i * xyzDistance[0] * 2);
                        result[1].Add(row * xyzDistance[1] * 2);
                        result[2].Add(layer * xyzDistance[2]);
                    }
                }
            }
        }

        return result;
    }
    
    public static List<List<int>> GenerateNt(List<List<double>> akt, List<int> xyzQuantity) // 20*nel глобальний номер і-ого локального вузла. 
    {
        var totalQuantity = xyzQuantity[0] * xyzQuantity[1] * xyzQuantity[2]; // nel кількість елементів, "квадратів"
        var edgeRowPointsCount = 1 + 2 * xyzQuantity[0]; // кількість точок в повному рядку 
        var middleRowPointsCount = 1 + xyzQuantity[0]; // кількість точок в "половинчастому" рядку
        var fullLayerPointsCount = ((edgeRowPointsCount + middleRowPointsCount) * xyzQuantity[1]) + edgeRowPointsCount; // кількість точок в повному шарі = рядок+(рядок+половинка)*y. тобто якщо у=1 => 2 повні рядки та половинка
        var smallLayerPointsCount = middleRowPointsCount * (xyzQuantity[1] + 1); // кількість точок в маленькому, "половинчастому" шарі
        var bothLayersPointsCount = fullLayerPointsCount + smallLayerPointsCount; // сума в 2 шарах
        var bothRowsPointsCount = edgeRowPointsCount + middleRowPointsCount; // сума в двох рядках

        var nt = new List<List<int>>();
        
        for (var i = 0; i < 20; i++)
        {
            nt.Add(new List<int>(Enumerable.Repeat(0, totalQuantity)));
        } // заповнюємо нулями 

        var cubesInLayer = xyzQuantity[0] * xyzQuantity[1]; // кількість кубів в шарі 

        for (var z = 0; z < xyzQuantity[2]; z++)
        {
            for (var y = 0; y < xyzQuantity[1]; y++)
            {
                for (var x = 0; x < xyzQuantity[0]; x++)
                {
                    var cubeNumber = cubesInLayer * z + xyzQuantity[0] * y + x; // номер куба = кількість кубів в шарі*z + кільість кубів в рядку * у + х
                    nt[0][cubeNumber] = z * bothLayersPointsCount + y * bothRowsPointsCount + x * 2; // z* 2 шари + y* 2 рядки + х*2
                    nt[1][cubeNumber] = nt[0][cubeNumber] + 2; // минулий + 2 
                    nt[2][cubeNumber] = nt[1][cubeNumber] + bothRowsPointsCount; // перший + кількість в 2 рядках
                    nt[3][cubeNumber] = nt[2][cubeNumber] - 2; //минулий - 2 

                    for (var i = 4; i < 8; i++)
                    {
                        nt[i][cubeNumber] = nt[i - 4][cubeNumber] + bothLayersPointsCount; // верхні точки, просто добавляємо 2 шари для кожної 
                    }
                    // нижні серединки
                    nt[8][cubeNumber] = nt[0][cubeNumber] + 1; // перший +1
                    nt[9][cubeNumber] = z * bothLayersPointsCount + y * bothRowsPointsCount + edgeRowPointsCount + x + 1; // z* 2 шари + y* 2 рядки +кількість в повному рядку +х + 1
                    nt[10][cubeNumber] = nt[8][cubeNumber] + bothRowsPointsCount; // перша нижня серединка + 2 шари
                    nt[11][cubeNumber] = nt[9][cubeNumber] - 1; // 8 - 1 
                    //середини в середньмоу шарі
                    nt[12][cubeNumber] = z * bothLayersPointsCount + fullLayerPointsCount + y * middleRowPointsCount + x; // z* 2 шари + повний шар + середні рядки + х
                    nt[13][cubeNumber] = nt[12][cubeNumber] + 1; // минулий +1
                    nt[14][cubeNumber] = nt[13][cubeNumber] + middleRowPointsCount; // минулий + кількість в середньому рядку
                    nt[15][cubeNumber] = nt[14][cubeNumber] - 1; // минулий -1 
                    // верхні серединки 
                    for (var i = 16; i < 20; i++)
                    {
                        nt[i][cubeNumber] = nt[i - 8][cubeNumber] + bothLayersPointsCount; // нижні серединки + 2 шари  
                    }
                }
            }
        }

        return nt;
    }


}
using System;
using System.Collections.Generic;

public class Constants
{
    public static readonly List<double> C = new() { 5.0 / 9, 8.0 / 9, 5.0 / 9 };
    public static readonly List<double> Xg = new() { -Math.Sqrt(0.6), 0, Math.Sqrt(0.6) };

    
    public static readonly IReadOnlyDictionary<int, List<double>> NodeNumberToLocalCoords = new Dictionary<int, List<double>>
    {
        { 0, new List<double> { -1, -1, -1 } },
        { 1, new List<double> { 1, -1, -1 } },
        { 2, new List<double> { 1, 1, -1 } },
        { 3, new List<double> { -1, 1, -1 } },
        { 4, new List<double> { -1, -1, 1 } },
        { 5, new List<double> { 1, -1, 1 } },
        { 6, new List<double> { 1, 1, 1 } },
        { 7, new List<double> { -1, 1, 1 } },
        { 8, new List<double> { 0, -1, -1 } },
        { 9, new List<double> { 1, 0, -1 } },
        { 10, new List<double> { 0, 1, -1 } },
        { 11, new List<double> { -1, 0, -1 } },
        { 12, new List<double> { -1, -1, 0 } },
        { 13, new List<double> { 1, -1, 0 } },
        { 14, new List<double> { 1, 1, 0 } },
        { 15, new List<double> { -1, 1, 0 } },
        { 16, new List<double> { 0, -1, 1 } },
        { 17, new List<double> { 1, 0, 1 } },
        { 18, new List<double> { 0, 1, 1 } },
        { 19, new List<double> { -1, 0, 1 } }
    };
    
    public static readonly IReadOnlyDictionary<int, List<double>> NodeNumberToLocalCoords2d = new Dictionary<int, List<double>>
    {
        { 0, new List<double> { -1, -1 } },
        { 1, new List<double> { 1, -1 } },
        { 2, new List<double> { 1, 1 } },
        { 3, new List<double> { -1, 1 } },
        { 4, new List<double> { 0, -1 } },
        { 5, new List<double> { 1, 0 } },
        { 6, new List<double> { 0, 1 } },
        { 7, new List<double> { -1, 0 } }
    };
    
    public static readonly IReadOnlyDictionary<int, List<int>> SideToPoint = new Dictionary<int, List<int>>
    {
        { 1, new List<int> { 0, 1, 5, 4, 8, 13, 16, 12 } },
        { 2, new List<int> { 1, 2, 6, 5, 9, 14, 17, 13 } },
        { 3, new List<int> { 2, 3, 7, 6, 10, 15, 18, 14 } },
        { 4, new List<int> { 3, 0, 4, 7, 11, 12, 19, 15 } },
        { 5, new List<int> { 3, 2, 1, 0, 10, 9, 8, 11 } },
        { 6, new List<int> { 4, 5, 6, 7, 16, 17, 18, 19 } }
    };

}
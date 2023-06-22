using System.Collections.Generic;

public static class GaussNodes
{
    public static List<List<double>> GetTwentySevenGaussianNodes()
    {
        var res = new List<List<double>>();

        for (var i = 0; i < 3; i++)
        {
            for (var j = 0; j < 3; j++)
            {
                for (var k = 0; k < 3; k++)
                {
                    res.Add(new List<double> { Constants.Xg[i], Constants.Xg[j], Constants.Xg[k] });
                }
            }
        }

        return res;
    }

    public static List<List<double>> GetNineGaussianNodes()
    {
        var res = new List<List<double>>();

        for (var i = 0; i < 3; i++)
        {
            for (var j = 0; j < 3; j++)
            {
                res.Add(new List<double> { Constants.Xg[i], Constants.Xg[j] });
            }
        }

        return res;
    }
    
}
using System;

public static class DeltaPhi
{
    public static double DeltaPhiAlpha_8(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.125 * (1 + beta * bi) * (1 + gamma * gi) * (ai * (2 * alpha * ai + beta * bi + gamma * gi - 1));
    }

    public static double DeltaPhiAlpha_12(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.25 * (1 + beta * bi) * (1 + gamma * gi) * (ai - 2 * alpha * Math.Pow(bi, 2) * Math.Pow(gi, 2) - 3 * Math.Pow(alpha, 2) * ai * Math.Pow(bi, 2) * Math.Pow(gi, 2) - ai * Math.Pow(beta * ai * gi, 2) - ai * Math.Pow(gamma * ai * bi, 2));
    }

    public static double DeltaPhiBeta_8(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.125 * (1 + alpha * ai) * (1 + gamma * gi) * (bi * (alpha * ai + 2 * beta * bi + gamma * gi - 1));
    }

    public static double DeltaPhiBeta_12(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.25 * (1 + alpha * ai) * (1 + gamma * gi) * (bi - bi * Math.Pow(alpha * bi * gi, 2) - 2 * beta * Math.Pow(ai, 2) * Math.Pow(gi, 2) - 3 * Math.Pow(beta, 2) * Math.Pow(ai, 2) * bi * Math.Pow(gi, 2) - bi * Math.Pow(gamma * ai * bi, 2));
    }

    public static double DeltaPhiGamma_8(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.125 * (1 + alpha * ai) * (1 + beta * bi) * (gi * (alpha * ai + beta * bi + 2 * gamma * gi - 1));
    }

    public static double DeltaPhiGamma_12(double alpha, double beta, double gamma, int nodeLocalNumber)
    {
        var coords = Constants.NodeNumberToLocalCoords[nodeLocalNumber];
        var ai = coords[0];
        var bi = coords[1];
        var gi = coords[2];
        return 0.25 * (1 + alpha * ai) * (1 + beta * bi) * (gi * (1 - Math.Pow(alpha * bi * gi, 2) - Math.Pow(beta * ai * gi, 2)) - 2 * gamma * Math.Pow(ai, 2) * Math.Pow(bi, 2) - 3 * Math.Pow(gamma, 2) * gi * Math.Pow(ai, 2) * Math.Pow(bi, 2));
    }
}
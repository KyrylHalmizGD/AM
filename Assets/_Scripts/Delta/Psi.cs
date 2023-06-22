using System;
using System.Collections;

public static class Psi
{
    private static double Psi4(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return 0.25 * (1 + eta * etai) * (1 + tau * taui) * (eta * etai + tau * taui - 1);
    }

    private static  double Psi57(double eta, double tau, int nodeLocalNumber)
    {
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return 0.5 * (1 - eta * eta) * (1 + tau * taui);
    }

    private static  double Psi68(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        return 0.5 * (1 - tau * tau) * (1 + eta * etai);
    }

    public static double Psii(double eta, double tau, int nodeLocalNumber, int i)
    {
        if (((IList)new[] { 0, 1, 2, 3 }).Contains(i))
            return Psi4(eta, tau, nodeLocalNumber);
        if (((IList)new[] { 4, 6 }).Contains(i))
            return Psi57(eta, tau, nodeLocalNumber);
        if (((IList)new[] { 5, 7 }).Contains(i))
            return Psi68(eta, tau, nodeLocalNumber);

        throw new ArgumentException("Invalid value of i.");
    }
}
using System;

public static class DeltaPsi
{
    public static double DeltaPsiEta4(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return 0.25 * (1 + tau * taui) * (tau * taui * etai + 2 * eta * Math.Pow(etai, 2));
    }

    public static double DeltaPsiTau4(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return 0.25 * (1 + eta * etai) * (eta * etai * taui + 2 * tau * Math.Pow(taui, 2));
    }

    public static double DeltaPsiEta57(double eta, double tau, int nodeLocalNumber)
    {
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return -eta * (1 + tau * taui);
    }

    public static double DeltaPsiTau57(double eta, double tau, int nodeLocalNumber)
    {
        double taui = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][1];
        return 0.5 * taui * (1 - Math.Pow(eta, 2));
    }

    public static double DeltaPsiEta68(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        return 0.5 * etai * (1 - Math.Pow(tau, 2));
    }

    public static double DeltaPsiTau68(double eta, double tau, int nodeLocalNumber)
    {
        double etai = Constants.NodeNumberToLocalCoords2d[nodeLocalNumber][0];
        return -tau * (1 + eta * etai);
    }
}
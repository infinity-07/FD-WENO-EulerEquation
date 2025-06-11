#ifndef __UTILITY__HPP__
#define __UTILITY__HPP__
#include "../third_party/scTools.h"

constexpr double C1_3 = 1.0 / 3.0;
constexpr double C1_4 = 1.0 / 4.0;
constexpr double C1_6 = 1.0 / 6.0;
constexpr double C5_6 = 5.0 / 6.0;
constexpr double C7_6 = 7.0 / 6.0;
constexpr double C11_6 = 11.0 / 6.0;
constexpr double C13_12 = 13.0 / 12.0;

inline double WENO5threconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num)
{
    // SHU, C.-W. (1997).
    // Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws.

    double tmp1(0), tmp2(0);
    double beta0(0), beta1(0), beta2(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0);
    double h0(0), h1(0), h2(0), ph(0);

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // use nonlinear weights
    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    beta0 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta1 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    beta2 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    const double eps(1.0e-6);
    w0 = d0 / ((eps + beta0) * (eps + beta0));
    w1 = d1 / ((eps + beta1) * (eps + beta1));
    w2 = d2 / ((eps + beta2) * (eps + beta2));

    double sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    switch (gp_num)
    {
    case 0:
        h0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        h1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        h2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        h0 = (-1.0 / 60.0 - sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave;
        h1 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5) / 20.0) * uavep;
        h2 = (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave + (1.0 / 30.0 - sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavepp;
        break;
    case 2:
        h0 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave;
        h1 = (-1.0 / 60 - sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavep;
        h2 = (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave + (1.0 / 30.0 + sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5) / 20.0) * uavepp;
        break;
    case 3:
        h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        h1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    ph = w0 * h0 + w1 * h1 + w2 * h2;

    return ph;
}
inline double WENO5Zthreconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num)
{
    // Borges, R., Carmona, M., Costa, B., & Don, W. S. (2008).
    // An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws.
    // Journal of Computational Physics, 227(6), 3191-3211. https://doi.org/10.1016/j.jcp.2007.11.038

    double tmp1(0), tmp2(0);
    double betam(0), beta(0), betap(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0), sum(0);
    double p0(0), p1(0), p2(0), ph(0);

    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    betam = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    betap = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // wenoz
    const double eps(1.0e-20);
    double tau = fabs(betam - betap);
    w0 = d0 * (1 + pow(tau / (betam + eps), 2));
    w1 = d1 * (1 + pow(tau / (beta + eps), 2));
    w2 = d2 * (1 + pow(tau / (betap + eps), 2));

    sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    switch (gp_num)
    {
    case 0:
        p0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        p1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        p2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        p0 = (-1.0 / 60.0 - sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave;
        p1 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave + (1.0 / 30.0 - sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavepp;
        break;
    case 2:
        p0 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave;
        p1 = (-1.0 / 60 - sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave + (1.0 / 30.0 + sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5) / 20.0) * uavepp;
        break;
    case 3:
        p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}
inline double WENO5ZPtheconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num, double deltaX)
{
    // Acker, F., B. de R. Borges, R., & Costa, B. (2016).
    // An improved WENO-Z scheme.
    // Journal of Computational Physics, 313, 726-753. https://doi.org/10.1016/j.jcp.2016.01.038

    double tmp1(0), tmp2(0);
    double betam(0), beta(0), betap(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0), sum(0);
    double p0(0), p1(0), p2(0), ph(0);

    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    betam = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    betap = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    double tau = fabs(betam - betap);
    double lambda = pow(deltaX, 2 / 3);

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // wenozp
    const double eps(1.0e-40);
    w0 = d0 * (1 + pow((tau + eps) / (betam + eps), 2) + lambda * (betam + eps) / (tau + eps));
    w1 = d1 * (1 + pow((tau + eps) / (beta + eps), 2) + lambda * (beta + eps) / (tau + eps));
    w2 = d2 * (1 + pow((tau + eps) / (betap + eps), 2) + lambda * (betap + eps) / (tau + eps));

    sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    switch (gp_num)
    {
    case 0:
        p0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        p1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        p2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        p0 = (-1.0 / 60.0 - sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave;
        p1 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave + (1.0 / 30.0 - sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavepp;
        break;
    case 2:
        p0 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave;
        p1 = (-1.0 / 60 - sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave + (1.0 / 30.0 + sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5) / 20.0) * uavepp;
        break;
    case 3:
        p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}
inline double WENO5ZPItheconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, int gp_num)
{
    // Luo, X., &Wu, S.- p.(2021).
    // Improvement of the WENO - Z + scheme.
    // Computers &Fluids, 218. https : // doi.org/10.1016/j.compfluid.2021.104855

    double tmp1(0), tmp2(0);
    double betam(0), beta(0), betap(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0), sum(0);
    double p0(0), p1(0), p2(0), ph(0);

    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    betam = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    betap = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    double tau = fabs(betam - betap);
    double lambda = 3.3;
    double betaMAX = std::max(betam, std::max(beta, betap));
    double betaMIN = std::min(betam, std::min(beta, betap));
    double alpha = 1 - (betaMIN / (betaMAX + 1e-40));

    // use linear weights
    switch (gp_num)
    {
    case 0:
        d0 = 0.3;
        d1 = 0.6;
        d2 = 0.1;
        break;
    case 1:
        d0 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 2:
        d0 = (91.0 + 9.0 * sqrt(5.0)) / 440.0;
        d1 = 129.0 / 220.0;
        d2 = (91.0 - 9.0 * sqrt(5.0)) / 440.0;
        break;
    case 3:
        d0 = 0.1;
        d1 = 0.6;
        d2 = 0.3;
        break;
    default:
        std::cout << gp_num << std::endl;
        std::cout << "check the gauss points for weno reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    // wenozpi
    const double eps(1.0e-40);
    w0 = d0 * (1 + pow(tau / (betam + eps), 2) + lambda * alpha * (betam) / (betaMAX + eps));
    w1 = d1 * (1 + pow(tau / (beta + eps), 2) + lambda * alpha * (beta) / (betaMAX + eps));
    w2 = d2 * (1 + pow(tau / (betap + eps), 2) + lambda * alpha * (betap) / (betaMAX + eps));

    sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    switch (gp_num)
    {
    case 0:
        p0 = -1.0 / 6.0 * uavemm + 5.0 / 6.0 * uavem + 1.0 / 3.0 * uave;
        p1 = 1.0 / 3.0 * uavem + 5.0 / 6.0 * uave - 1.0 / 6 * uavep;
        p2 = 11.0 / 6.0 * uave - 7.0 / 6.0 * uavep + 1.0 / 3.0 * uavepp;
        break;
    case 1:
        p0 = (-1.0 / 60.0 - sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 + sqrt(5) / 5.0) * uavem + (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave;
        p1 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60 - sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave + (1.0 / 30.0 - sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavepp;
        break;
    case 2:
        p0 = (-1.0 / 60.0 + sqrt(5) / 20.0) * uavemm + (1.0 / 30.0 - sqrt(5) / 5.0) * uavem + (59.0 / 60.0 + 3 * sqrt(5) / 20) * uave;
        p1 = (-1.0 / 60 - sqrt(5) / 20.0) * uavem + (31.0 / 30.0) * uave + (-1.0 / 60.0 + sqrt(5) / 20.0) * uavep;
        p2 = (59.0 / 60.0 - 3 * sqrt(5) / 20.0) * uave + (1.0 / 30.0 + sqrt(5) / 5.0) * uavep + (-1.0 / 60.0 - sqrt(5) / 20.0) * uavepp;
        break;
    case 3:
        p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
        p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
        p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;
        break;
    default:
        std::cout << "Check the Gauss points for WENO reconstruction!" << std::endl;
        std::cin.get();
        exit(1);
    }

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}

inline double WENO5threconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp)
{
    // SHU, C.-W. (1997).
    // Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws.

    // use linear weights
    const double d0 = 0.1;
    const double d1 = 0.6;
    const double d2 = 0.3;

    // use nonlinear weights
    // smooth indicator for current component
    double tmp1 = (uavemm - 2 * uavem + uave);
    double tmp2 = (uavemm - 4 * uavem + 3 * uave);
    const double beta0 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    const double beta1 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    const double beta2 = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    const double eps(1.0e-6);
    double w0 = d0 / ((eps + beta0) * (eps + beta0));
    double w1 = d1 / ((eps + beta1) * (eps + beta1));
    double w2 = d2 / ((eps + beta2) * (eps + beta2));

    const double sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    double h0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
    double h1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
    double h2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;

    return w0 * h0 + w1 * h1 + w2 * h2;
}
inline double WENO5Zthreconstruction(double &uavemm, double &uavem, double &uave, double &uavep, double &uavepp)
{
    // Borges, R., Carmona, M., Costa, B., & Don, W. S. (2008).
    // An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws.
    // Journal of Computational Physics, 227(6), 3191-3211. https://doi.org/10.1016/j.jcp.2007.11.038

    // smooth indicator for each stencil
    double tmp1(0), tmp2(0);
    tmp1 = (uavemm - 2.0 * uavem + uave);
    tmp2 = (uavemm - 4.0 * uavem + 3.0 * uave);
    const double betam = C13_12 * tmp1 * tmp1 + C1_4 * tmp2 * tmp2;

    tmp1 = (uavem - 2.0 * uave + uavep);
    tmp2 = (uavem - uavep);
    const double beta = C13_12 * tmp1 * tmp1 + C1_4 * tmp2 * tmp2;

    tmp1 = (uave - 2.0 * uavep + uavepp);
    tmp2 = (3.0 * uave - 4.0 * uavep + uavepp);
    const double betap = C13_12 * tmp1 * tmp1 + C1_4 * tmp2 * tmp2;

    const double eps(1.0e-20); // wenoz
    const double inv_betam = 1.0 / (betam + eps);
    const double inv_beta = 1.0 / (beta + eps);
    const double inv_betap = 1.0 / (betap + eps);

    // use linear weights
    const double d0 = 0.1;
    const double d1 = 0.6;
    const double d2 = 0.3;

    // wenoz nonlinear weights
    const double tau = fabs(betam - betap);
    double w0 = d0 * (1.0 + (tau * inv_betam * tau * inv_betam));
    double w1 = d1 * (1.0 + (tau * inv_beta * tau * inv_beta));
    double w2 = d2 * (1.0 + (tau * inv_betap * tau * inv_betap));

    const double weightSum = w0 + w1 + w2;
    w0 /= weightSum;
    w1 /= weightSum;
    w2 /= weightSum;

    // p(x)
    const double p0 = C1_3 * uavemm - C7_6 * uavem + C11_6 * uave;
    const double p1 = -C1_6 * uavem + C5_6 * uave + C1_3 * uavep;
    const double p2 = C1_3 * uave + C5_6 * uavep - C1_6 * uavepp;

    return w0 * p0 + w1 * p1 + w2 * p2;
}
inline double WENO5ZPtheconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp, double deltaX)
{
    // Acker, F., B. de R. Borges, R., & Costa, B. (2016).
    // An improved WENO-Z scheme.
    // Journal of Computational Physics, 313, 726-753. https://doi.org/10.1016/j.jcp.2016.01.038

    double tmp1(0), tmp2(0);
    double betam(0), beta(0), betap(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0), sum(0);
    double p0(0), p1(0), p2(0), ph(0);

    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    betam = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    betap = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    double tau = fabs(betam - betap);
    double lambda = pow(deltaX, 2 / 3);

    // use linear weights
    d0 = 0.1;
    d1 = 0.6;
    d2 = 0.3;

    // wenozp
    const double eps(1.0e-40);
    w0 = d0 * (1 + pow((tau + eps) / (betam + eps), 2) + lambda * (betam + eps) / (tau + eps));
    w1 = d1 * (1 + pow((tau + eps) / (beta + eps), 2) + lambda * (beta + eps) / (tau + eps));
    w2 = d2 * (1 + pow((tau + eps) / (betap + eps), 2) + lambda * (betap + eps) / (tau + eps));

    sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
    p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
    p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}
inline double WENO5ZPItheconstruction(double uavemm, double uavem, double uave, double uavep, double uavepp)
{
    // Luo, X., &Wu, S.- p.(2021).
    // Improvement of the WENO - Z + scheme.
    // Computers &Fluids, 218. https : // doi.org/10.1016/j.compfluid.2021.104855

    double tmp1(0), tmp2(0);
    double betam(0), beta(0), betap(0);
    double d0(0), d1(0), d2(0);
    double w0(0), w1(0), w2(0), sum(0);
    double p0(0), p1(0), p2(0), ph(0);

    // smooth indicator for current component
    tmp1 = (uavemm - 2 * uavem + uave);
    tmp2 = (uavemm - 4 * uavem + 3 * uave);
    betam = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uavem - 2 * uave + uavep);
    tmp2 = (uavem - uavep);
    beta = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    tmp1 = (uave - 2 * uavep + uavepp);
    tmp2 = (3 * uave - 4 * uavep + uavepp);
    betap = 13.0 / 12.0 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2;

    double tau = fabs(betam - betap);
    double lambda = 3.3;
    double betaMAX = std::max(betam, std::max(beta, betap));
    double betaMIN = std::min(betam, std::min(beta, betap));
    double alpha = 1 - (betaMIN / (betaMAX + 1e-40));

    // use linear weights
    d0 = 0.1;
    d1 = 0.6;
    d2 = 0.3;

    // wenozpi
    const double eps(1.0e-40);
    w0 = d0 * (1 + pow(tau / (betam + eps), 2) + lambda * alpha * (betam) / (betaMAX + eps));
    w1 = d1 * (1 + pow(tau / (beta + eps), 2) + lambda * alpha * (beta) / (betaMAX + eps));
    w2 = d2 * (1 + pow(tau / (betap + eps), 2) + lambda * alpha * (betap) / (betaMAX + eps));

    sum = w0 + w1 + w2;

    w0 = w0 / sum;
    w1 = w1 / sum;
    w2 = w2 / sum;

    // p(x)
    p0 = 1.0 / 3.0 * uavemm - 7.0 / 6.0 * uavem + 11.0 / 6.0 * uave;
    p1 = -1.0 / 6 * uavem + 5.0 / 6.0 * uave + 1.0 / 3.0 * uavep;
    p2 = 1.0 / 3.0 * uave + 5.0 / 6.0 * uavep - 1.0 / 6.0 * uavepp;

    ph = w0 * p0 + w1 * p1 + w2 * p2;

    return ph;
}

#endif
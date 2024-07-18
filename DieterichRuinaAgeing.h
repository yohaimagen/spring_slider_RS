#ifndef DIETERICHRUINAAGEING_20201027_H
#define DIETERICHRUINAAGEING_20201027_H

#include "Zero.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iostream>

#include <math.h> 


class DieterichRuinaAgeing {
public:


    double V0;
    double b;
    double f0;
    double a;
    double eta;
    double L;
    double sn;
    double Vinit;
    double Vp;
    double k;
    double yield_point_init;



    

    double psi_init(double tau) const {
        if (Vinit == 0.0) {
            return f0;
        }
        double s = sinh((tau - eta * Vinit) / (a * sn));
        double l = log((2.0 * V0 / Vinit) * s);
        return a * l;
    }


    auto slip_rate(double tau, double psi) const -> double{
        double V = 0.0;
        double alpha = 0.0;
        double beta = tau / eta;
        if (alpha > beta)
        {
            std::swap(alpha, beta);
        }
        auto fF = [this, &tau, &psi](double V) -> double {
            double sn = this->sn;
            double V0 = this->V0;
            double a = this->a;
            return tau - this->F(sn, V, psi, a, V0) - eta * V;
        };
        try {
            V = zeroIn(alpha, beta, fF);
        } catch (std::exception const&) {
            std::cout << "sigma_n = " << sn << std::endl
                        << "|tau| = " << tau << std::endl
                        << "psi = " << psi << std::endl
                        << "L = " << a << std::endl
                        << "U = " << b << std::endl
                        << "F(L) = " << fF(a) << std::endl
                        << "F(U) = " << fF(b) << std::endl;
            throw;

        }
        return V;
    }

    double state_rhs(double V, double psi) const {
        return b * V0 / L * (exp((f0 - psi) / b) - V / V0);
    }

   


//private:

    double F(double sn, double V, double psi, double a, double V0) const {
        double e = exp(psi / a);
        double f = a * asinh((V / (2.0 * V0)) * e);
        double tau = sn * f;
        return tau;
    }

    double Finv(double snAbs, double tauAbs, double psi, double a, double V0) const {
        // We have
        // V = 2 V_0 sinh(tauAbs / (a snAbs)) exp(-psi / a)
        // substituting r = tauAbs / snAbs and sinh(x) = (exp(x) - exp(-x)) / 2 we obtain
        // V = V_0 (exp((r - psi) / a) - exp(-(r + psi) / a))
        double r = tauAbs / snAbs;
        return V0 * (exp((r - psi) / a) - exp(-(r + psi) / a));
    }

};

#endif // DIETERICHRUINAAGEING_20201027_H

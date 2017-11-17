/*************************************************************************
 *
 *  This file is part of HFHSL.
 *
 *  Author: Ivo Filot <i.a.w.filot@tue.nl>
 *
 *  HFHSL is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  HFHSL is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HFHSL.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#include <Eigen/Eigenvalues>

#include "cgf.h"
#include "integrals.h"

/*
 * Calculates the electronic energy for a H atom using the STO-3G and STO-6G basis sets.
 *
 * All calculations are performed in standard units.
 */

int main() {

    /*
     * STO-3G
     */

    // construct the STO-3G CGF for H
    const vec3 pos1(0.0, 0.0, 0.0);
    CGF cgf1(pos1);
    cgf1.add_gto(CGF::GTO_S, 3.4252509099999999,  0.15432897000000001, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.62391373000000006, 0.53532813999999995, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.16885539999999999, 0.44463454000000002, pos1);

    // create integrator object
    Integrator integrator;

    // calculate integrals
    double S = integrator.overlap(cgf1, cgf1);
    double T = integrator.kinetic(cgf1, cgf1);
    double V = integrator.nuclear(cgf1, cgf1, pos1, 1.0);

    std::cout << "Overlap integral:    " << S << std::endl;
    std::cout << "Kinetic integral:    " << T << std::endl;
    std::cout << "Nuclear integral:    " << V << std::endl;
    std::cout << "Total energy:        " << T + V << std::endl;

    std::cout << "-----------------------------" << std::endl;

    /*
     * STO-6G
     */

    // construct the STO-6G CGF for H
    CGF cgf2(pos1);
    cgf2.add_gto(CGF::GTO_S, 35.523221, 0.009164, pos1);
    cgf2.add_gto(CGF::GTO_S, 6.513144, 0.049361, pos1);
    cgf2.add_gto(CGF::GTO_S, 1.822143, 0.168538, pos1);
    cgf2.add_gto(CGF::GTO_S, 0.625955, 0.370563, pos1);
    cgf2.add_gto(CGF::GTO_S, 0.243077, 0.416492, pos1);
    cgf2.add_gto(CGF::GTO_S, 0.100112, 0.130334, pos1);

    // calculate integrals
    S = integrator.overlap(cgf2, cgf2);
    T = integrator.kinetic(cgf2, cgf2);
    V = integrator.nuclear(cgf2, cgf2, pos1, 1.0);

    std::cout << "Overlap integral:    " << S << std::endl;
    std::cout << "Kinetic integral:    " << T << std::endl;
    std::cout << "Nuclear integral:    " << V << std::endl;
    std::cout << "Total energy:        " << T + V << std::endl;

    std::cout << "-----------------------------" << std::endl;

    return 0;
}

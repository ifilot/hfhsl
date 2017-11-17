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

#include "cgf.h"
#include "integrals.h"

/*
 * Calculate the electronic energy for the He atom
 *
 * All calculations are performed in standard units.
 */

int main() {

    // construct the CGF for He
    const vec3 pos1(0.0, 0.0, 0.0);
    CGF cgf1(pos1);
    cgf1.add_gto(CGF::GTO_S, 6.3624213899999997,  0.15432897000000001, pos1);
    cgf1.add_gto(CGF::GTO_S, 1.1589229999999999, 0.53532813999999995, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.31364978999999998, 0.44463454000000002, pos1);

    // create integrator object
    Integrator integrator;

    // calculate the overlap, kinetic and nuclear integrals
    double S = integrator.overlap(cgf1, cgf1);
    double T = integrator.kinetic(cgf1, cgf1);
    double V = integrator.nuclear(cgf1, cgf1, pos1, 2.0);   // NOTE THE 2.0 HERE! This is He, not H.

    // calculate the two-electron repulsion integral
    double te = integrator.repulsion(cgf1, cgf1, cgf1, cgf1);

    // evaluate the electronic energy of He
    double E = 2 * (T + V) + te;

    std::cout << "Overlap integral:               " << S << std::endl;
    std::cout << "Kinetic integral:               " << T << std::endl;
    std::cout << "Nuclear integral:               " << V << std::endl;
    std::cout << "Electron repulsion integral:    " << te << std::endl;
    std::cout << "Total energy:                   " << E << std::endl;

    return 0;
}

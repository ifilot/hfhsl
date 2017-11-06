/*************************************************************************
 *
 *  This file is part of HFCXX.
 *
 *  Author: Ivo Filot <i.a.w.filot@tue.nl>
 *
 *  HFCXX is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  HFCXX is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HFCXX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

#include <Eigen/Eigenvalues>

#include "cgf.h"
#include "integrals.h"

/*
 * Calculates the energy of H2 using the Hartree-Fock Self-Consistent Field
 * method. An STO-3G basis set is used in the description of the 1s orbitals
 * of the two H atoms. The H atoms are positioned 1.4 a.u. apart.
 *
 * All calculations are performed in standard units.
 */

int main() {

    // construct two CGFs for two H atoms 1.4 a.u. apart
    const vec3 pos1(0.0, 0.0, 0.0);
    CGF cgf1(pos1);
    cgf1.add_gto(CGF::GTO_S, 6.3624213899999997,  0.15432897000000001, pos1);
    cgf1.add_gto(CGF::GTO_S, 1.1589229999999999, 0.53532813999999995, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.31364978999999998, 0.44463454000000002, pos1);

    // Collect all CGFS in a vector object
    std::vector<CGF> cgfs;
    cgfs.push_back(cgf1);

    // create integrator object
    Integrator integrator;

    // calculate the integral values using the integrator class
    double S = integrator.overlap(cgf1, cgf1);
    double T = integrator.kinetic(cgf1, cgf1);
    double V = integrator.nuclear(cgf1, cgf1, pos1, 2.0);

    // calculate all two-electron integrals
    double te = integrator.repulsion(cgf1, cgf1, cgf1, cgf1);

    // Calculate 1-electron Hamiltonian matrix
    double E = 2 * (T + V) + te;

    std::cout << "Overlap integral:               " << S << std::endl;
    std::cout << "Kinetic integral:               " << T << std::endl;
    std::cout << "Nuclear integral:               " << V << std::endl;
    std::cout << "Electron repulsion integral:    " << te << std::endl;
    std::cout << "Total energy:                   " << E << std::endl;

    return 0;
}

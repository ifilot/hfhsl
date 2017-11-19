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
    cgf1.add_gto(CGF::GTO_S, 3.4252509099999999,  0.15432897000000001, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.62391373000000006, 0.53532813999999995, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.16885539999999999, 0.44463454000000002, pos1);

    const vec3 pos2(0.0, 0.0, 1.4);
    CGF cgf2(pos2);
    cgf2.add_gto(CGF::GTO_S, 3.4252509099999999,  0.15432897000000001, pos2);
    cgf2.add_gto(CGF::GTO_S, 0.62391373000000006, 0.53532813999999995, pos2);
    cgf2.add_gto(CGF::GTO_S, 0.16885539999999999, 0.44463454000000002, pos2);

    // Collect all CGFS in a vector object
    std::vector<CGF> cgfs;
    cgfs.push_back(cgf1);
    cgfs.push_back(cgf2);

    // create integrator object
    Integrator integrator;

    // Construct 2x2 matrices to hold values for the overlap,
    // kinetic and two nuclear integral values, respectively.
    Eigen::Matrix2d S, T, V1, V2;

    // calculate the integral values using the integrator class
    for(unsigned int i=0; i<2; i++) {
        for(unsigned int j=0; j<2; j++) {
            S(i,j) = integrator.overlap(cgfs[i], cgfs[j]);
        }
    }

    // print overlap matrix
    std::cout << "S-Matrix" << std::endl;
    std::cout << S << std::endl;
    std::cout << std::endl;

    // perform matrix diagonalization
    Eigen::EigenSolver<Eigen::Matrix2d> es(S, true);
    Eigen::Matrix2d D = es.eigenvalues().real().asDiagonal();
    Eigen::Matrix2d U = es.eigenvectors().real();

    // Output diagonal matrix
    std::cout << "D-Matrix" << std::endl;
    std::cout << D << std::endl;
    std::cout << std::endl;

    // Output unitary matrix
    std::cout << "U-Matrix" << std::endl;
    std::cout << U << std::endl;
    std::cout << std::endl;

    return 0;
}

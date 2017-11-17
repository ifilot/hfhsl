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
            T(i,j) = integrator.kinetic(cgfs[i], cgfs[j]);
            V1(i,j) = integrator.nuclear(cgfs[i], cgfs[j], pos1, 1.0);
            V2(i,j) = integrator.nuclear(cgfs[i], cgfs[j], pos2, 1.0);
        }
    }

    // calculate all two-electron integrals
    unsigned int size = integrator.teindex(2,2,2,2);
    std::vector<double> tedouble(size, -1.0);
    for(unsigned int i=0; i<2; i++) {
        for(unsigned int j=0; j<2; j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<2; k++) {
                for(unsigned int l=0; l<2; l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = integrator.teindex(i,j,k,l);

                        // this avoids recalculating an integral which has
                        // already been evaluated
                        if(tedouble[idx] != -1.0) {
                            continue;
                        }

                        tedouble[idx] = integrator.repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l]);
                    }
                }
            }
        }
    }

    // print overlap matrix
    std::cout << "S-Matrix" << std::endl;
    std::cout << S << std::endl;
    std::cout << std::endl;

    // print kinetic matrix
    std::cout << "T-Matrix" << std::endl;
    std::cout << T << std::endl;
    std::cout << std::endl;

    // print nuclear attracton matrix
    std::cout << "V1-Matrix" << std::endl;
    std::cout << V1 << std::endl;
    std::cout << std::endl;

    // print nuclear attracton matrix
    std::cout << "V2-Matrix" << std::endl;
    std::cout << V2 << std::endl;
    std::cout << std::endl;

    // print all unique two electorn integrals
    std::cout << "Two electron integrals:" << std::endl;
    for(unsigned int i=0; i<2; i++) {
        for(unsigned int j=0; j<2; j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<2; k++) {
                for(unsigned int l=0; l<2; l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = integrator.teindex(i,j,k,l);
                        std::cout << "[" << idx << "] " << "(" << i+1 << "," << j+1 << "," << k+1 << "," << l+1 << ") ";
                        std::cout << tedouble[idx] << std::endl;
                    }
                }
            }
        }
    }

    return 0;
}

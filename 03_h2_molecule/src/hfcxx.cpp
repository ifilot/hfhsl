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
                        tedouble[idx] = integrator.repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l]);
                    }
                }
            }
        }
    }

    // perform a canonical diagonalization on the overlap matrix to
    // obtain orthonormal spinorbitals (required for the Slater
    // Determinant)
    Eigen::EigenSolver<Eigen::Matrix2d> es(S, true);
    Eigen::Matrix2d D;
    D(0,0) = 1.0 / std::sqrt(es.eigenvalues()(0,0).real() );
    D(1,1) = 1.0 / std::sqrt(es.eigenvalues()(1,0).real() );

    Eigen::Matrix2d U = es.eigenvectors().real();

    // Calculate the transformation matrix
    Eigen::Matrix2d X = U * D;
    Eigen::Matrix2d Xp = X.transpose();

    // Create Density matrices
    Eigen::Matrix2d P;
    Eigen::Matrix2d Pnew;

    // Calculate 1-electron Hamiltonian matrix
    Eigen::Matrix2d H = T + V1 + V2;

    // Create two-electron Hamiltonian matrix
    Eigen::Matrix2d G;

    // define mixing parameter in the SCF iterations
    static const double alpha = 0.5;

    // hold energy value of previous iteration
    double energy_old = 0.0;

    // difference between previous and current SCF iteration
    double energy_difference = 1.0;

    // keep track of number of iterations
    unsigned int loop_counter = 0;

    /*
     * START ITERATIVE PROCEDURE
     */
    while(energy_difference > 1e-5) {
        loop_counter++; // increment loop counter

        // Populate two-electron hamiltonian matrix
        for(unsigned int i=0; i<2; i++) {
            for(unsigned int j=0; j<2; j++) {
                G(i,j) = 0.; /* reset G matrix */
                for(unsigned int k=0; k<2; k++) {
                    for(unsigned int l=0; l<2; l++) {
                        unsigned int index1 = integrator.teindex(i,j,l,k);
                        unsigned int index2 = integrator.teindex(i,k,l,j);
                        G(i,j) += P(k,l) * (tedouble[index1] - 0.5 * tedouble[index2]);
                    }
                }
            }
        }

        // Calculate Fock Matrix
        Eigen::Matrix2d F = H + G;

        // Transform Fock Matrix using our basis transformation matrix
        Eigen::Matrix2d Fp = Xp * F * X;

        // Calculate eigenvalues and vectors
        es.compute(Fp, true);
        Eigen::Matrix2d Cc = es.eigenvectors().real();

        // Calculate energy
        double energy = 0.0;
        Eigen::Matrix2d M = H + F;
        for(unsigned int i=0; i<2; i++) {
            for(unsigned int j=0; j<2; j++) {
                energy += P(j,i) * M(i,j);
            }
        }

        // Add nuclear repulsion to the orbital energies
        energy = energy * 0.5 + 1.0 / 1.4;

        // Obtain true coefficient matrix using the transformation matrix
        Eigen::Matrix2d C = X * Cc;

        Pnew = Eigen::Matrix2d::Zero();
        for(unsigned int i=0; i<2; i++) {
            for(unsigned int j=0; j<2; j++) {
                for(unsigned int k=0; k<1; k++) {
                    Pnew(i,j) += 2.0 * C(i,k) * C(j,k);
                }
            }
        }

        for(unsigned int i=0; i<2; i++) {
            for(unsigned int j=0; j<2; j++) {
                P(i,j) = (1.0-alpha) * Pnew(i,j) + alpha * P(i,j);
            }
        }

        // report energy and calculate difference with previous energy value
        energy_difference = std::abs(energy - energy_old);
        energy_old = energy;

        std::cout << loop_counter << ":\t" << energy << std::endl;
    }

    std::cout << "Stopping because energy convergence was achieved." << std::endl;

    return 0;
}

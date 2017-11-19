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
#include "sort.h"

/*
 * Calculates the energy of H2 using the Hartree-Fock Self-Consistent Field
 * method. An STO-3G basis set is used in the description of the 1s orbitals
 * of the two H atoms. The H atoms are positioned 1.4 a.u. apart.
 *
 * All calculations are performed in standard units.
 */

int main() {

    /*********************************
     *
     * STEP 1: Define nuclei and basis functions
     *
     *********************************/

    std::vector<CGF> cgfs;
    // construct the CGFs for CO
    // C
    const vec3 pos1(0.0, 0.0, 0.0);
    // 1S
    cgfs.emplace_back(pos1);
    cgfs.back().add_gto(CGF::GTO_S, 71.616837,  0.154329, pos1);
    cgfs.back().add_gto(CGF::GTO_S, 13.045096, 0.535328, pos1);
    cgfs.back().add_gto(CGF::GTO_S, 3.530512, 0.444635, pos1);
    // 2S
    cgfs.emplace_back(pos1);
    cgfs.back().add_gto(CGF::GTO_S, 2.941249,  -0.099967, pos1);
    cgfs.back().add_gto(CGF::GTO_S, 0.683483, 0.399513, pos1);
    cgfs.back().add_gto(CGF::GTO_S, 0.222290, 0.700115, pos1);
    // 2Px
    cgfs.emplace_back(pos1);
    cgfs.back().add_gto(CGF::GTO_PX, 2.941249,  0.155916, pos1);
    cgfs.back().add_gto(CGF::GTO_PX, 0.683483, 0.607684, pos1);
    cgfs.back().add_gto(CGF::GTO_PX, 0.222290, 0.391957, pos1);
    // 2Py
    cgfs.emplace_back(pos1);
    cgfs.back().add_gto(CGF::GTO_PY, 2.941249,  0.155916, pos1);
    cgfs.back().add_gto(CGF::GTO_PY, 0.683483, 0.607684, pos1);
    cgfs.back().add_gto(CGF::GTO_PY, 0.222290, 0.391957, pos1);
    // 2Pz
    cgfs.emplace_back(pos1);
    cgfs.back().add_gto(CGF::GTO_PZ, 2.941249,  0.155916, pos1);
    cgfs.back().add_gto(CGF::GTO_PZ, 0.683483, 0.607684, pos1);
    cgfs.back().add_gto(CGF::GTO_PZ, 0.222290, 0.391957, pos1);

    // O
    const vec3 pos2(0.0, 0.0, 2.116);
    // 1S
    cgfs.emplace_back(pos2);
    cgfs.back().add_gto(CGF::GTO_S, 130.709320,  0.154329, pos2);
    cgfs.back().add_gto(CGF::GTO_S, 23.808861, 0.535328, pos2);
    cgfs.back().add_gto(CGF::GTO_S, 6.443608, 0.444635, pos2);
    // 2S
    cgfs.emplace_back(pos2);
    cgfs.back().add_gto(CGF::GTO_S, 5.033151,  -0.099967, pos2);
    cgfs.back().add_gto(CGF::GTO_S, 1.169596, 0.399513, pos2);
    cgfs.back().add_gto(CGF::GTO_S, 0.380389, 0.700115, pos2);
    // 2Px
    cgfs.emplace_back(pos2);
    cgfs.back().add_gto(CGF::GTO_PX, 5.033151,  0.155916, pos2);
    cgfs.back().add_gto(CGF::GTO_PX, 1.169596, 0.607684, pos2);
    cgfs.back().add_gto(CGF::GTO_PX, 0.380389, 0.391957, pos2);
    // 2Py
    cgfs.emplace_back(pos2);
    cgfs.back().add_gto(CGF::GTO_PY, 5.033151,  0.155916, pos2);
    cgfs.back().add_gto(CGF::GTO_PY, 1.169596, 0.607684, pos2);
    cgfs.back().add_gto(CGF::GTO_PY, 0.380389, 0.391957, pos2);
    // 2Pz
    cgfs.emplace_back(pos2);
    cgfs.back().add_gto(CGF::GTO_PZ, 5.033151,  0.155916, pos2);
    cgfs.back().add_gto(CGF::GTO_PZ, 1.169596, 0.607684, pos2);
    cgfs.back().add_gto(CGF::GTO_PZ, 0.380389, 0.391957, pos2);

    // create integrator object
    Integrator integrator;

    /*********************************
     *
     * STEP 2: Calculate S, T, V, H and TE-integrals
     *
     *********************************/

    // Construct 2x2 matrices to hold values for the overlap,
    // kinetic and two nuclear integral values, respectively.
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());
    Eigen::MatrixXd V1 = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());
    Eigen::MatrixXd V2 = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());

    // calculate the integral values using the integrator class
    for(unsigned int i=0; i<cgfs.size(); i++) {
        for(unsigned int j=0; j<cgfs.size(); j++) {
            S(i,j) = integrator.overlap(cgfs[i], cgfs[j]);
            T(i,j) = integrator.kinetic(cgfs[i], cgfs[j]);
            V1(i,j) = integrator.nuclear(cgfs[i], cgfs[j], pos1, 6.0);
            V2(i,j) = integrator.nuclear(cgfs[i], cgfs[j], pos2, 8.0);
        }
    }

    // Calculate 1-electron Hamiltonian matrix
    Eigen::MatrixXd H = T + V1 + V2;

    // calculate all two-electron integrals
    unsigned int size = integrator.teindex(cgfs.size(),cgfs.size(),cgfs.size(),cgfs.size());
    std::vector<double> tedouble(size, -1.0);
    for(unsigned int i=0; i<cgfs.size(); i++) {
        for(unsigned int j=0; j<cgfs.size(); j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<cgfs.size(); k++) {
                for(unsigned int l=0; l<cgfs.size(); l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = integrator.teindex(i,j,k,l);
                        tedouble[idx] = integrator.repulsion(cgfs[i], cgfs[j], cgfs[k], cgfs[l]);
                    }
                }
            }
        }
    }

    /*********************************
     *
     * STEP 3: Calculate transformation matrix
     *
     *********************************/

    // perform a canonical diagonalization on the overlap matrix to
    // obtain orthonormal spinorbitals (required for the Slater
    // Determinant)
    Eigen::EigenSolver<Eigen::MatrixXd> es(S, true);
    Eigen::MatrixXd D = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd U = es.eigenvectors().real();
    sort_eigenvalues(U,D);
    for(unsigned int i=0; i<cgfs.size(); i++) {
        D(i,i) = 1.0 / sqrt(D(i,i));
    }

    std::cout << S << std::endl;

    // Calculate the transformation matrix
    Eigen::MatrixXd X = U * D;
    Eigen::MatrixXd Xp = X.transpose();

    /*********************************
     *
     * STEP 4: Obtain initial guess for density matrix
     *
     *********************************/

    // Create Density matrices
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());
    Eigen::MatrixXd Pnew = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());

    // Create two-electron Hamiltonian matrix
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());

    // define mixing parameter in the SCF iterations
    static const double alpha = 0.5;

    // hold energy value of previous iteration
    double energy_old = 0.0;

    // difference between previous and current SCF iteration 
    // (initialize with some large number)
    double energy_difference = 1.0;

    // keep track of number of iterations
    unsigned int loop_counter = 0;

    // construct object to store eigenvectors and eigenvalues of the Fock-matrix in
    Eigen::MatrixXd C;
    Eigen::VectorXd orbital_energies;

    /*
     * START ITERATIVE PROCEDURE
     */
    while(energy_difference > 1e-5 && loop_counter < 100) {
        loop_counter++; // increment loop counter

        /*********************************
         *
         * STEP 5: Calculate G, H, F and F' from P
         *
         *********************************/

        // Populate two-electron hamiltonian matrix
        for(unsigned int i=0; i<cgfs.size(); i++) {
            for(unsigned int j=0; j<cgfs.size(); j++) {
                G(i,j) = 0.; /* reset G matrix */
                for(unsigned int k=0; k<cgfs.size(); k++) {
                    for(unsigned int l=0; l<cgfs.size(); l++) {
                        unsigned int index1 = integrator.teindex(i,j,l,k);
                        unsigned int index2 = integrator.teindex(i,k,l,j);
                        G(i,j) += P(k,l) * (tedouble[index1] - 0.5 * tedouble[index2]);
                    }
                }
            }
        }

        // Calculate Fock Matrix
        Eigen::MatrixXd F = H + G;

        // Transform Fock Matrix using our basis transformation matrix
        Eigen::MatrixXd Fp = Xp * F * X;

        /*********************************
         *
         * STEP 6: Diagonalize F' to obtain C' and e
         *
         *********************************/

        // Calculate eigenvalues and vectors
        es.compute(Fp, true);
        Eigen::MatrixXd Cc = es.eigenvectors().real();
        Eigen::MatrixXd en = es.eigenvalues().real().asDiagonal();

        // Calculate energy
        double energy = 0.0;
        Eigen::MatrixXd M = H + F;
        for(unsigned int i=0; i<cgfs.size(); i++) {
            for(unsigned int j=0; j<cgfs.size(); j++) {
                energy += P(j,i) * M(i,j);
            }
        }

        // Add nuclear repulsion to the orbital energies
        energy = energy * 0.5 + 6.0 * 8.0 / pos2[2];

        /*********************************
         *
         * STEP 7: Diagonalize C from C'
         *
         *********************************/

        // Obtain true coefficient matrix using the transformation matrix
        C = X * Cc;
        sort_eigenvalues(C,en);
        orbital_energies = en.diagonal();

        /*********************************
         *
         * STEP 8: Calculate new P from C
         *
         *********************************/

        // obtain new P matrix from the old matrix
        Pnew = Eigen::MatrixXd::Zero(cgfs.size(), cgfs.size());
        for(unsigned int i=0; i<cgfs.size(); i++) {
            for(unsigned int j=0; j<cgfs.size(); j++) {
                for(unsigned int k=0; k<7; k++) {
                    Pnew(i,j) += 2.0 * C(i,k) * C(j,k);
                }
            }
        }

        // construct new P-matrix by mixing old and new matrix. This is not always necessary, 
        // but sometimes the SCF calculation does not converge without it.
        for(unsigned int i=0; i<cgfs.size(); i++) {
            for(unsigned int j=0; j<cgfs.size(); j++) {
                P(i,j) = (1.0-alpha) * Pnew(i,j) + alpha * P(i,j);
            }
        }

        // report energy and calculate difference with previous energy value
        energy_difference = std::abs(energy - energy_old);
        energy_old = energy;

        std::cout << loop_counter << ":\t" << energy << std::endl;
    }

    std::cout << "Stopping because energy convergence was achieved." << std::endl;

    // outputting orbital energies
    for(unsigned int i=0; i<10; i++) {
        std::cout << "e" << i+1 << "= " << orbital_energies[i] << std::endl;
    }

    return 0;
}

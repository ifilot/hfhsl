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

#include "density_plotter.h"

/**
 * @brief      default constructor
 */
DensityPlotter::DensityPlotter() {}

/**
 * @brief      export electron densities of each MO to a CHGCAR file (VASP density format)
 *
 * @param[in]  cgfs  reference to vector of contracted gaussian functionals
 * @param[in]  C     reference to coefficient matrix C
 */
void DensityPlotter::plot_densities_chgcar(const std::vector<CGF>& cgfs, const Eigen::MatrixXd& C) {
	static const double range = 10;
	static const double inc = 0.2;
	double volume = range * range * range * 8;

	for(unsigned int i=0; i<cgfs.size(); i++) {
		Eigen::VectorXd coefs = C.col(i);

		std::ofstream outfile((boost::format("%04i.dat") % (i+1)).str());

		outfile << "DENSITY MO" << std::endl;
		outfile << "1.0" << std::endl;
		outfile << "20 0 0" << std::endl;
		outfile << "0 20 0" << std::endl;
		outfile << "0 0 20" << std::endl;
		outfile << "1" << std::endl;
		outfile << "Direct" << std::endl;
		outfile << "0 0 0" << std::endl;
		outfile << "" << std::endl;
		outfile << (int)(range/inc*2+1) << " ";
		outfile << (int)(range/inc*2+1) << " ";
		outfile << (int)(range/inc*2+1) << std::endl;

		unsigned int cnt = 0;
		for(double x=-range; x<=range; x+=inc) {
			for(double y=-range; y<=range; y+=inc) {
				for(double z=-range; z<=range; z+=inc) {

					double sum = 0.0;
					for(unsigned int j=0; j<cgfs.size(); j++) {
						sum += coefs[j] * cgfs[j].get_value(x,y,z);
					}

					outfile << (sum * sum) * 2.0 * volume << " ";

					cnt++;
					if(cnt % 5 == 0) {
						outfile << std::endl;
						cnt = 0;
					}
				}
			}
		}
	}
}

/**
 * @brief      export electron densities of each MO to a cube file (Gaussian
 *             format)
 *
 * @param[in]  atoms  vector of atom number of atoms
 * @param[in]  pos    position of the atoms
 * @param[in]  cgfs   set of contracted gaussian functions (basis set)
 * @param[in]  C      coefficient matrix
 */
void DensityPlotter::plot_densities_cube(const std::vector<unsigned int>& atoms, 
										 const std::vector<vec3>& pos,
										 const std::vector<CGF>& cgfs, 
										 const Eigen::MatrixXd& C) {
	static const double range = 10;
	static const double inc = 0.2;

	for(unsigned int i=0; i<cgfs.size(); i++) {
		Eigen::VectorXd coefs = C.col(i);

		std::string filename = (boost::format("%04i.cube") % (i+1)).str();
		std::cout << "Outputting cube file to: " << filename << std::endl;
		std::ofstream outfile(filename);

		outfile << "HFHSL CUBE FILE." << std::endl;
 		outfile << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << std::endl;
    	// print number of atoms and center of unit cell
    	outfile << (boost::format("%i  %g  %g  %g") % atoms.size() % range % range % range).str() << std::endl;

    	// print unit cell information
    	outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % inc % 0.0 % 0.0).str() << std::endl;
    	outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % 0.0 % inc % 0.0).str() << std::endl;
    	outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % 0.0 % 0.0 % inc).str() << std::endl;

    	// print atom positions
    	for(unsigned int i=0; i<atoms.size(); i++) {
    		outfile << (boost::format("%i  %g  %g  %g") % atoms[i] % pos[i][0] % pos[i][1] % pos[i][2]).str() << std::endl;
    	}

		unsigned int cnt = 0;
		for(double x=-range; x<=range; x+=inc) {
			for(double y=-range; y<=range; y+=inc) {
				for(double z=-range; z<=range; z+=inc) {

					double sum = 0.0;
					for(unsigned int j=0; j<cgfs.size(); j++) {
						sum += coefs[j] * cgfs[j].get_value(x,y,z);
					}

					outfile << (boost::format("  %g") % ((sum * sum) * 2.0)).str();

					cnt++;
					if(cnt % 6 == 0) {
						outfile << std::endl;
						cnt = 0;
					}
				}
				outfile << std::endl;
				cnt = 0;
			}
		}
	}
}
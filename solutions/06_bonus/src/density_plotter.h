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

#ifndef _DENSITY_PLOTTER_H
#define _DENSITY_PLOTTER_H

#include <fstream>
#include <string>
#include <Eigen/Eigenvalues>
#include <boost/format.hpp>

#include "gamma.h"
#include "cgf.h"

class DensityPlotter {
private:

public:
	DensityPlotter();

	void plot_densities_chgcar(const std::vector<CGF>& cgfs, const Eigen::MatrixXd& C);

	void plot_densities_cube(const std::vector<unsigned int>& atoms, 
							 const std::vector<vec3>& pos,
							 const std::vector<CGF>& cgfs, 
							 const Eigen::MatrixXd& C);

private:

};

#endif // _DENSITY_PLOTTER_H
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

#ifndef _HFHSL_SORT_H
#define _HFHSL_SORT_H

#include <Eigen/Eigenvalues>

#include <vector>
#include <utility>
#include <algorithm>

/**
 * @brief      sort the eigenvalues and eigenvectors
 *
 * @param      U     Unitary eigenvector matrix
 * @param      D     Diagonal eigenvalue matrix
 */
void sort_eigenvalues(Eigen::MatrixXd& U, Eigen::MatrixXd& D) {
    std::vector<std::pair<double, Eigen::VectorXd>> e;
    for(unsigned int i=0; i<U.cols(); i++) {
        if(U(0,i) > 0.0) {
            e.push_back(std::pair<double, Eigen::VectorXd>(D(i,i), -U.col(i)));
        } else{
            e.push_back(std::pair<double, Eigen::VectorXd>(D(i,i), U.col(i)));
        }
    }
    std::sort(e.begin(), e.end(), [](auto &left, auto &right) {
        return left.first < right.first;
    });

    unsigned int j=0;
    for(const auto& p : e) {
        for(unsigned int i=0; i<U.cols(); i++) {
            U(i,j) = p.second(i);
        }
        D(j,j) = p.first;
        j++;
    }
}

#endif // _HFHSL_SORT_H
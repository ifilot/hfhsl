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

#ifndef _CGF_H
#define _CGF_H

#include <iostream>
#include <Eigen/Dense>
#include <boost/math/special_functions/factorials.hpp>

typedef Eigen::Vector3d vec3;

/*
 * Gaussian Type Orbital
 *
 * N * (x-X)^l * (y-Y)^m * (z-Z)^n * exp(-alpha * r^2)
 *
 * where r = sqrt(x^2 + y^2 + z^2)
 * and N a normalization constant such that <GTO|GTO> = 1
 *
 */

class GTO { // Gaussian Type Orbital
private:
    double c;               // coefficient
    double alpha;           // alpha value in the exponent
    unsigned int l,m,n;     // powers of the polynomial
    vec3 r;                 // position vector (unit = Bohr)
    double norm;            // normalization constant

public:
    GTO(double _c,
        const vec3& _r,
        double _alpha,
        unsigned int _l,
        unsigned int _m,
        unsigned int _n);

    inline const double get_alpha() const {
        return this->alpha;
    }

    inline const unsigned int get_l() const {
        return this->l;
    }

    inline const unsigned int get_m() const {
        return this->m;
    }

    inline const unsigned int get_n() const {
        return this->n;
    }

    inline const double get_norm() const {
        return this->norm;
    }

    inline const double get_coefficient() const {
        return this->c;
    }

    inline const vec3& get_r() const {
        return this->r;
    }

private:
    void calculate_normalization_constant();
};


class CGF { // Contracted Gaussian Function
private:
    std::vector<GTO> gtos;  // vector holding all gtos
    vec3 r;                 // position of the CGF

public:
    CGF(const vec3& _r);

    enum{
        GTO_S,
        GTO_PX,
        GTO_PY,
        GTO_PZ,
        GTO_DX2,
        GTO_DXY,
        GTO_DXZ,
        GTO_DY2,
        GTO_DYZ,
        GTO_DZ2,

        NUM_GTO
    };

    inline const unsigned int size() const {
        return this->gtos.size();
    }

    inline const double get_norm_gto(const unsigned int i) const {
        return this->gtos[i].get_norm();
    }

    inline const double get_coefficient_gto(const unsigned int i) const {
        return this->gtos[i].get_coefficient();
    }

    inline const GTO& get_gto(const unsigned int i) const {
        return this->gtos[i];
    }

    // add a Gaussian type orbital to the CGF
    void add_gto(unsigned int type,  // type of the orbital (see above for the list)
                 double alpha,       // alpha value
                 double c,           // coefficient
                 const vec3& vec3);  // position
};

#endif //_CGF_H

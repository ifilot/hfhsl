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

#include "cgf.h"

GTO::GTO(double _c,
         const vec3& _r,     // position (unit = Bohr)
         double _alpha,
         unsigned int _l,
         unsigned int _m,
         unsigned int _n):
    c(_c),
    alpha(_alpha),
    l(_l),
    m(_m),
    n(_n),
    r(_r) {

    // calculate the normalization constant
    this->calculate_normalization_constant();
}

void GTO::calculate_normalization_constant() {
    static const double pi = 3.14159265359;

    double nom =   std::pow(2.0, 2.0 * (l + m + n) + 3.0 / 2.0) *
                   std::pow(alpha, (l + m + n) + 3.0 / 2.0);

    double denom = (l < 1 ? 1 : boost::math::double_factorial<double>(2 * l - 1) )*
                   (m < 1 ? 1 : boost::math::double_factorial<double>(2 * m - 1) )*
                   (n < 1 ? 1 : boost::math::double_factorial<double>(2 * n - 1) )*
                   std::pow(pi, 3.0 / 2.0);

    this->norm = std::sqrt(nom / denom);
}

double GTO::get_value(double x, double y, double z) const {
    return this->c * 
           std::pow(x - this->r[0], l) *
           std::pow(y - this->r[1], m) *
           std::pow(z - this->r[2], n) *
           std::exp(-this->alpha * (this->r - vec3(x,y,z)).squaredNorm());
}

CGF::CGF(const vec3& _r):
    r(_r) {
        // do nothing
}

double CGF::get_value(double x, double y, double z) const {
    double sum = 0.0;

    for(const auto& gto : this->gtos) {
        sum += gto.get_value(x,y,z);
    }

    return sum;
}

void CGF::add_gto(unsigned int type,  // type of the orbital (see above for the list)
                  double alpha,       // alpha value
                  double c,           // coefficient
                  const vec3& vec3) { // position

    switch(type) {
        // S ORBITAL
        case GTO_S:
            gtos.push_back(GTO(c, r, alpha, 0,0,0));
        break;

        // P ORBITALS
        case GTO_PX:
            gtos.push_back(GTO(c, r, alpha, 1,0,0));
        break;
        case GTO_PY:
            gtos.push_back(GTO(c, r, alpha, 0,1,0));
        break;
        case GTO_PZ:
            gtos.push_back(GTO(c, r, alpha, 0,0,1));
        break;

        // D ORBITALS
        case GTO_DX2:
            gtos.push_back(GTO(c, r, alpha, 2,0,0));
        break;
        case GTO_DXY:
            gtos.push_back(GTO(c, r, alpha, 1,1,0));
        break;
        case GTO_DXZ:
            gtos.push_back(GTO(c, r, alpha, 1,0,1));
        break;
        case GTO_DY2:
            gtos.push_back(GTO(c, r, alpha, 0,2,0));
        break;
        case GTO_DYZ:
            gtos.push_back(GTO(c, r, alpha, 0,1,1));
        break;
        case GTO_DZ2:
            gtos.push_back(GTO(c, r, alpha, 0,0,2));
        break;

        default:
            std::cerr << "Undefined orbital type. Exiting..." << std::endl;
            exit(-1);
        break;
    }
}

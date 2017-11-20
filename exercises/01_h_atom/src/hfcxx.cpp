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

#include "cgf.h"
#include "integrals.h"

/*
 * Calculates the electronic energy for a H atom using the STO-3G and STO-6G basis sets.
 *
 * All calculations are performed in standard units.
 */

int main() {

    // construct the STO-3G CGF for H
    const vec3 pos1(0.0, 0.0, 0.0); // position of the atom
    CGF cgf1(pos1);
    cgf1.add_gto(CGF::GTO_S, 3.4252509099999999,  0.15432897000000001, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.62391373000000006, 0.53532813999999995, pos1);
    cgf1.add_gto(CGF::GTO_S, 0.16885539999999999, 0.44463454000000002, pos1);

    // create integrator object
    Integrator integrator;

    /*
     * START CODING HERE
     */

    std::cout << "Done" << std::endl;

    return 0;
}

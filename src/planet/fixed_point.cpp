/*****************************************************************************
 *   Copyright (C) 2004-2015 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include "fixed_point.h"
#include "../epoch.h"
#include "../exceptions.h"

namespace kep_toolbox{ namespace planet {

fixedp::fixedp(double r1 = 0.0, 
               double r2 = 0.0, 
               double r3 = 0.0, 
               double v1 = 0.0, 
               double v2 = 0.0, 
               double v3 = 0.0) : base()
{
	fixed_r[0] = r1;   fixed_r[1] = r2;   fixed_r[2] = r3;
	fixed_v[0] = v1;   fixed_v[1] = v2;   fixed_v[2] = v3;
}

/// Polymorphic copy constructor.
planet_ptr fixedp::clone() const
{
	return planet_ptr(new fixedp(*this));
}

/// Returns the stored position and velocity
/// (considered constant for any given time)
void fixedp::eph_impl(double mjd2000, array3D &r, array3D &v) const {

	for (unsigned int i=0; i<3; i++) {
		r[i] = fixed_r[i];   v[i] = fixed_v[i];
	}
}

/// Extra informations streamed in human readable format
std::string fixedp::human_readable_extra() const {
	std::ostringstream s;
	s << "Ephemerides type: fixed coordinates" << std::endl;
	return s.str();
}

}} //namespace

// Serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(kep_toolbox::planet::fixedp)
// Serialization code (END)

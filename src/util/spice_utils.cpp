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

#include <string>
#include <sstream>

#include "spice_utils.h"

namespace kep_toolbox { namespace util {

/// Load SPICE kernel
/**
 * This function wraps the SPICE toolbox furnsh_c routine. You can find the original documentation
 * at http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
 *
 * You must call this function to load the corresponding kernel before calling the ephemerdies method of a
 * kep_toolbox::planet_spice, otherwise an exception will be thrown by the spice library.
 *
 * \param[in] file_name an std::string containing the filename of the kernel to be loaded (e.g. "de432s.bsp")
 *
 * \throws See original documentation for the possible exceptions
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

void load_spice_kernel(std::string file_name) {
	/// Transferring error handling from spice to kep_toolbox
	erract_c("SET",0,const_cast<char*>("RETURN") ) ;
	/// Loading the kernel
    furnsh_c ( file_name.c_str() );
    /// Handling errors
    if (failed_c()) {
    	std::ostringstream msg;
    	msg << "Kernel file " << file_name << " could not be loaded by SPICE" << std::endl;
        reset_c();
    	throw_value_error(msg.str());
    }
}

/// Transforms kep_toolbox epoch to SPICE epoch
/**
 * This function wraps the SPICE toolbox furnsh_c routine. You can find the original documentation
 * at http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
 *
 * You must call this function to load the corresponding kernel before calling the ephemerdies method of a
 * kep_toolbox::planet_spice, otherwise an exception will be thrown by the spice library.
 *
 * \param[in] ep a kep_toolbox epoch
 * \returns a SpiceDouble containing the SPICE epoch (seconds from 2000.01.01 12:00:00)
 *
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */
SpiceDouble epoch_to_spice(kep_toolbox::epoch ep) {
	return ( (ep.mjd2000() -.5) * ASTRO_DAY2SEC );
}

SpiceDouble epoch_to_spice(double mjd2000) {
    return ( (mjd2000 -.5) * ASTRO_DAY2SEC );
}

}} // namespaces


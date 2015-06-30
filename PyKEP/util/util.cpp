/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 3 of the License, or       *
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

// Workaround for http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#ifdef _WIN32
#include <cmath>
#endif

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>

#ifdef PYKEP_USING_SPICE
	#include "../../src/util/spice_utils.h"
#endif

using namespace boost::python;

BOOST_PYTHON_MODULE(_util) {
	// Disable docstring c++ signature to allow sphinx autodoc to work properly
	docstring_options doc_options;
	doc_options.disable_signatures();
	
#ifdef PYKEP_USING_SPICE
	// Spice utilities
	def("load_spice_kernel",&kep_toolbox::util::load_spice_kernel,
			  "PyKEP.util.load_spice_kernel(file_name)\n\n"
			  "- file_name: string containing the kernel file to load\n\n"
			  "Loads the SPICE kernel specified by the filename into memory. \n\n"
			  ".. note::\n\n"
			  "   The kernel must be in memory before its used, for example, when computing a planets.spice ephemerides\n\n"
			  "Example:: \n\n"
			  "  util.load_spice_kernel('de432s.bsp')"
	);
	def("inspect_kernels",&kep_toolbox::util::inspect_kernels,
			  "PyKEP.util.inspect_kernels()\n\n"
			  "Lists the SPICE kernels in memory that have been loaded with load_spice_kernel().\n\n"
	);
	def("unload_spice_kernel",&kep_toolbox::util::unload_spice_kernel,
			  "PyKEP.util.unload_spice_kernel(file_name)\n\n"
              "- file_name: string containing the kernel file to unload\n\n"
              "Unloads the SPICE kernel specified by the filename from memory. \n\n"
              ".. note::\n\n"
              "   It is safe to unload a kernel that hasn't loaded - in that case, no action is taken.\n\n"
			  "Example:: \n\n"
			  "  util.unload_spice_kernel('de432s.bsp')"
	);
	def("get_eph_coverage",&kep_toolbox::util::get_eph_coverage,
			  "PyKEP.util.get_eph_coverage(lsk, spk)\n\n"
			  "Reports the timespan for the ephemerides kernel.\n\n"
	);
#endif
}


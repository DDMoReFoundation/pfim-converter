/* ----------------------------------------------------------------------------
 * This file is part of PharmML to PFIM 4.0 Converter.  
 * 
 * Copyright (C) 2016 jointly by the following organisations:- 
 * 1. INSERM, Paris, France
 * 2. Universite Paris Diderot, Paris, France
 * 3. EMBL European Bioinformatics Institute (EBML-EBI), Hinxton, UK
 * 4. Cyprotex Discovery Ltd, Macclesfield, England, UK
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation. A copy of the license agreement is provided
 * in the file named "LICENSE.txt" included with this software distribution.
 * ----------------------------------------------------------------------------
 */

package inserm.converters.pfim.scriptlets;

import crx.converter.engine.scriptlets.Factory;
import crx.converter.engine.scriptlets.Host;

/**
 * Factory class to instantiate the scripting host handle.<br/>
 * The scripting syntax is python.
 */
public class FactoryImpl implements Factory {
	@Override
	public Host newInstance() {
		try {
			return new HostImpl();
		} catch (Exception e) {
			return null;
		}
	}
}

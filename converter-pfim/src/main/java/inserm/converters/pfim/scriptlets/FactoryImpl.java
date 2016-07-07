/*******************************************************************************
 * Copyright (C) 2015 Cyprotex Discovery Ltd - All rights reserved.
 ******************************************************************************/

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

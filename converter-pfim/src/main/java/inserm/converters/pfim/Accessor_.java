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

package inserm.converters.pfim;

import java.util.List;

import crx.converter.engine.Accessor;
import eu.ddmore.libpharmml.dom.PharmML;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.dataset.MapType;
import eu.ddmore.libpharmml.dom.dataset.TargetMapping;
import eu.ddmore.libpharmml.dom.trialdesign.DosingVariable;

/**
 * PFIM specific modification of the standard
 */
public class Accessor_ extends Accessor {
	/**
	 * Constructor
	 * @param dom_ Model Handle
	 */
	public Accessor_(PharmML dom_) {
		super(dom_);
	}

	/**
	 * Fetch the model element declared as a dosing variable.
	 * @param target
	 * @return PharmMLRootType
	 */
	public PharmMLRootType fetchElement(DosingVariable target) {
		if (target == null) return null;
		else if (target.getTargetMapping() != null) {
			TargetMapping tm = target.getTargetMapping();
			List<MapType> mappings = tm.getListOfMap();
			if (mappings.isEmpty()) throw new IllegalStateException("Target mapping list is empty, no dose variable specified.");
			else if (mappings.size() > 1) throw new IllegalStateException("Multiple dose targets not yet suppported for a single PharmML administration element.");
			else if (mappings.size() == 1) return fetchElement(mappings.get(0).getModelSymbol());
			else return null;
		}
		else if (target.getSymbRef() != null) return fetchElement(target.getSymbRef());
		else return null;
	}
}

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

package inserm.converters.pfim.tree;

import inserm.converters.pfim.CategoricalCovariateRef;
import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.tree.BaseTreeMaker;
import crx.converter.tree.BinaryTree;

/**
 * PFIM Tree Maker.
 */
public class TreeMaker_ extends BaseTreeMaker {
	/**
	 * Constructor
	 */
	public TreeMaker_() { super(); }
	
	private BinaryTree createTree(CategoricalCovariateRef cref) { return createRootTree(cref, "Categorical Covariate Reference"); }
	
	private BinaryTree createTree(IndividualParameterAssignment ipa) { return createRootTree(ipa, ""); }
	
	@Override
	public BinaryTree newInstance(Object o) {
		BinaryTree bt = null;
		flushNestedTreeReferences();
		
		if (o instanceof IndividualParameterAssignment) bt = createTree((IndividualParameterAssignment) o);
		else if (o instanceof CategoricalCovariateRef) bt = createTree((CategoricalCovariateRef) o);
		else bt = super.newInstance(o);
		
		return bt;
	}
}
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

package inserm.converters.pfim.parts;

import java.util.ArrayList;
import java.util.List;

import crx.converter.engine.Part;
import crx.converter.spi.ILexer;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;

/**
 * A part represents a block of code to be created from elements in a PharmMl model.
 *
 */
public abstract class PartImpl implements Part {
	private static Rhs rhs(Double value) {
		Rhs rhs = new Rhs();
		rhs.setScalar(new RealValue(value));
		return rhs;
	}
	
	/**
	 * Array of model elements assigned a matrices.
	 */
	protected List<Object> has_matrices = new ArrayList<Object>();
	
	/**
	 * The converter instance.
	 */
	protected ILexer lexer = null;
	
	/**
	 * A list of model symbols to be used in the model block.
	 */
	protected ArrayList<String> symbolIds = new ArrayList<String>();
	
	/**
	 * Flag that the syntax trees have been built.
	 */
	protected boolean treeBuilt = false;
	
	/**
	 * Method to build the syntax trees bound to a part.
	 */
	abstract public void buildTrees();
	
	/**
	 * Check the assignment statement of a parameter and set to zero if undefined.
	 * @param p Parameter
	 */
	protected void checkAssignment(PopulationParameter p) {
		if (p == null) return;
		boolean isAssignmentOkay = false;
		if (p.getAssign() != null) {
			Rhs rhs = p.getAssign();
			if (rhs.getContent() != null) isAssignmentOkay = true;
		}
		
		if (!isAssignmentOkay) p.setAssign(rhs(0.0));
	}
	
	/**
	 * Get the name of the code block
	 * @return java.lang.String
	 */
	abstract public String getName();
	
	/**
	 * Get a list of symbols bound to a part.
	 * @return java.util.List<java.lang.String>
	 */
	abstract public List<String> getSymbolIds();
	
	/**
	 * Check if the Part has a bound model symbol.
	 * @param name Name of a variable
	 * @return boolean
	 */
	abstract public boolean hasSymbolId(String name);
	
	/**
	 * Check if the model symbol is okay, i.e. not null and unique.<br/>
	 * The Lexer instance can be set to permit duplicate variables.
	 * @param symbolId
	 * @return boolean
	 * @see crx.converter.spi.ILexer#setDuplicateVariablesPermitted
	 */
	protected boolean isSymbolIdOkay(String symbolId) {
		boolean isOkay = true;
		if (symbolId == null) throw new IllegalStateException("Symbol ID cannot be NULL.");
		
		if (symbolIds.contains(symbolId) && !lexer.isDuplicateVariablesPermitted()) isOkay = false;
		else symbolIds.add(symbolId);

		return isOkay;
	}
}

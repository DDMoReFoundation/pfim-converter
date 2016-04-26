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

import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;

/**
 * A block of individual parameters.
 *
 */
public abstract class IndividualParameterBlockImpl extends PartImpl {
	private List<IndividualParameter> external_indiv_params = new ArrayList<IndividualParameter>();
	
	/**
	 * Individual parameters bound to a block.
	 */
	protected List<IndividualParameter> indiv_params = new ArrayList<IndividualParameter>();
	
	/**
	 * Add an individual parameters to the block
	 * @param ip Individual parameter
	 * @return boolean
	 */
	public boolean addIndividualParameter(IndividualParameter ip) {
		if (ip != null) {
			if (!indiv_params.contains(ip)) {
				indiv_params.add(ip);
				external_indiv_params.add(ip); // Just to keep track of IP from outside the current block.
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Build the parameter trees associated with the code.
	 */
	protected void buildIndividualParameterTrees() {
		TreeMaker tm = lexer.getTreeMaker();
		if (tm == null) throw new NullPointerException("Treemaker is NULL");
		
		for (IndividualParameter ip : indiv_params) {
			if (ip == null) continue;
			
			BinaryTree bt = tm.newInstance(ip);
			IndividualParameterAssignment assignment = new IndividualParameterAssignment(ip);
			bt.nodes.get(0).data = assignment;
			lexer.addStatement(ip, bt);
			lexer.updateNestedTrees();
		}
	}
	
	/**
	 * The the list of individual parameters.
	 * @return java.util.List<IndividualParameter>
	 */
	public List<IndividualParameter> getIndividualParameters() {
		return indiv_params;
	}
	
	/**
	 * Flag whether the individual parameter block is empty.
	 * @return boolean
	 */
	public boolean hasIndividualParameters() { return indiv_params.size() > 0; }
	
	/**
	 * Flag whether an parameter is an external accessed element, shared between structural models.
	 * @param ip
	 * @return boolean
	 */
	protected boolean isExternalIndividualParameter(IndividualParameter ip) {
		if (ip != null) return external_indiv_params.contains(ip);
		return false;
	}
}

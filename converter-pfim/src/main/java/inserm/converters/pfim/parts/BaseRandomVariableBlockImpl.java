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

import static crx.converter.engine.PharmMLTypeChecker.isRandomVariable;

import java.util.ArrayList;
import java.util.List;

import crx.converter.engine.common.CorrelationRef;
import crx.converter.spi.blocks.BaseRandomVariableBlock;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.NestedTreeRef;
import crx.converter.tree.Node;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.modeldefn.CorrelatedRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.Correlation;
import eu.ddmore.libpharmml.dom.modeldefn.Pairwise;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;

/**
 * Base code block that contains random variables.
 */
public abstract class BaseRandomVariableBlockImpl extends IndividualParameterBlockImpl implements BaseRandomVariableBlock  {
	/**
	 * A list of correlation references.
	 */
	protected List<CorrelationRef> correlations = new ArrayList<CorrelationRef>();
	
	/**
	 * Flag if a random variable is linked to a pairwise statement.
	 */
	protected List<ParameterRandomVariable> linked_rvs = new ArrayList<ParameterRandomVariable>();
	
	/**
	 * List of random variables declared in the block.
	 */
	protected List<ParameterRandomVariable> rvs = new ArrayList<ParameterRandomVariable>();
	
	/**
	 * Build the AST linked to a correlation reference
	 * @param corr Correlation Reference
	 */
	protected void buildTrees(CorrelationRef corr) {
		if (corr == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		
		BinaryTree bt_ = new BinaryTree();
		Node root_node = new Node(corr);
		root_node.root = true;
		
		lexer.addStatement(corr, bt_);	 
		if (corr.is_pairwise) {
			if (corr.rnd1 != null) {
				lexer.addStatement(corr.rnd1, tm.newInstance(corr.rnd1));
				lexer.updateNestedTrees();
			}
			
			if (corr.rnd2 != null) {
				lexer.addStatement(corr.rnd2, tm.newInstance(corr.rnd2));
				lexer.updateNestedTrees();
			}
			
			// Make note of parameter references in the statements.
			List<BinaryTree> trees = new ArrayList<BinaryTree>();
			if (corr.correlationCoefficient != null) {
				BinaryTree bt = tm.newInstance(corr.correlationCoefficient);
				trees.add(bt);
				for (NestedTreeRef ntr : tm.getNestedTrees()) {
					if (ntr == null) continue;
					trees.add(ntr.bt);
				}
				lexer.addStatement(corr.correlationCoefficient, bt);
				lexer.updateNestedTrees();
			}
			
			if (corr.covariance != null) {
				BinaryTree bt = tm.newInstance(corr.covariance);
				trees.add(bt);
				for (NestedTreeRef ntr : tm.getNestedTrees()) {
					if (ntr == null) continue;
					trees.add(ntr.bt);
				}
				lexer.addStatement(corr.covariance, bt);	
				lexer.updateNestedTrees();
			}
			
			corr.findParameterReferences(trees);
		}
	}
	
	/**
	 * Get a list of correlation references.
	 * @return java.util.List<CorrelationRef>
	 */
	public List<CorrelationRef> getCorrelations() {
		return correlations;
	}
	
	/**
	 * A list of pairwise linked variables
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable>
	 */
	public List<ParameterRandomVariable> getLinkedRandomVariables() { return linked_rvs; }
	
	/**
	 * Create the random variable terms linked to a PharmML correlation.
	 * @param corr Correlation term
	 */
	protected void register(Correlation corr) {
		if (corr == null) return;
		if (corr.getMatrix() != null) {
			CorrelationRef corr_ref = new CorrelationRef();
			corr_ref.matrix = corr.getMatrix(); 
		} else if (corr.getPairwise() != null) {
			Pairwise pw = corr.getPairwise();
			
			CorrelatedRandomVariable rnd1 = pw.getRandomVariable1();
			CorrelatedRandomVariable rnd2 = pw.getRandomVariable2();
			
			if (rnd1 != null && rnd2 != null) {
				PharmMLRootType v1 = null, v2 = null;
				
				if (rnd1.getSymbRef() != null) v1 = lexer.getAccessor().fetchElement(rnd1.getSymbRef());
				if (rnd2.getSymbRef() != null) v2 = lexer.getAccessor().fetchElement(rnd2.getSymbRef());
				
				if (isRandomVariable(v1) && isRandomVariable(v2)) {
					ParameterRandomVariable rv1 = (ParameterRandomVariable) v1;
					ParameterRandomVariable rv2 = (ParameterRandomVariable) v2;
					
					if (rvs.contains(rv1)) { 
						rvs.remove(rv1);
						linked_rvs.add(rv1);
					}
					if (rvs.contains(rv2)) {
						rvs.remove(rv2);
						linked_rvs.add(rv2);
					}
					
					correlations.add(new CorrelationRef(lexer, rv1, rv2, pw.getCorrelationCoefficient(), pw.getCovariance())); 
				} else { 
					throw new IllegalStateException("Random variable reference of wrong type.");
				}
			} else 
				throw new NullPointerException("Correlation random variable is NULL.");
		}
	}
}

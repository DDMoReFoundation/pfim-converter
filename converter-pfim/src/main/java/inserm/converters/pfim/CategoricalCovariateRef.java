/*******************************************************************************
 * Copyright (C) 2016 Cyprotex Discovery Ltd - All rights reserved.
 ******************************************************************************/

package inserm.converters.pfim;

import java.util.ArrayList;
import java.util.List;

import crx.converter.spi.ILexer;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.modeldefn.CategoricalCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.Category;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;

/**
 * Place the categorical covariate data into a structure as opposed to the containment hierachy as defined within PharmML.
 * Mostly used to specify and construct a declaration term in generated code.
 */
public class CategoricalCovariateRef {
	public String blkId = null;
	private ILexer c = null;
	public List<String> categories = new ArrayList<String>();
	public CovariateDefinition cov = null;
	
	/**
	 * Constructor
	 * @param lexer_ Converter
	 * @param cov_ Covariate Declaration
	 * @param blkId_ Bulk Identifier of the Covariate Model
	 */
	public CategoricalCovariateRef(ILexer lexer_, CovariateDefinition cov_, String blkId_) {
		if (lexer_ == null) throw new NullPointerException("The converter reference cannot be NULL.");
		if (cov_ == null) throw new NullPointerException("A categorical covariate cannot be NULL.");
		if (blkId_ == null) throw new NullPointerException("The covariate model block identifier cannot be NULL.");
		cov = cov_;
		c = lexer_;
	}
	
	/**
	 * Build the AST required by the category definitions.
	 */
	public void buildTrees() {
		CategoricalCovariate catcov = cov.getCategorical();
		TreeMaker tm = c.getTreeMaker();
		
		for (Category category : catcov.getListOfCategory()) {
			if (category == null) continue;
			if (category.getCatId() == null) throw new IllegalStateException("A category ID is NULL. (cov='" + cov.getSymbId() + "')");
			
			c.addStatement(category, tm.newInstance(category));
			
			categories.add(category.getCatId());
			Rhs prob = category.getProbability();
			if (prob != null) {
				c.addStatement(prob, tm.newInstance(prob));
				c.updateNestedTrees();
			}
		}	
	}
	
	/**
	 * Checks if the reference contains the specified covariate.
	 * @param c
	 * @return boolean
	 */
	public boolean contains(CovariateDefinition c) {
		if (c == null) return false;
		else return cov.equals(c);
	}
	
	/**
	 * Get the name of the categorical covariate.
	 * @return String
	 */
	public String getName() { return cov.getSymbId(); }
}

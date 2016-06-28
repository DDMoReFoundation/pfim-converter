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

import static crx.converter.engine.PharmMLTypeChecker.isCommonParameter;
import static crx.converter.engine.PharmMLTypeChecker.isCovariate;
import inserm.converters.pfim.CategoricalCovariateRef;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.engine.common.CovariateParameterRef;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.CovariateBlock;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.modeldefn.CommonParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ContinuousCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateModel;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateTransformation;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;

/**
 * List of covariates read from the PharmML model.<br/>
 */
public class CovariateBlockImpl extends PartImpl implements CovariateBlock {
	private boolean categorical = false;
	private	Map<String, CategoricalCovariateRef> categorical_cov_map = new HashMap<String, CategoricalCovariateRef>();
	private List<String> categorical_cov_names = new ArrayList<String>();
	private List<CategoricalCovariateRef> categorical_covariates = new ArrayList<CategoricalCovariateRef>();
	private CovariateModel cm = null;
	private List<ContinuousCovariate> continuous_covariates = new ArrayList<ContinuousCovariate>();
	private Map<CovariateDefinition, CovariateParameterRef> cov_param_map = new HashMap<CovariateDefinition, CovariateParameterRef>();
	private List<CovariateDefinition> covariates = new ArrayList<CovariateDefinition>();
	private List<PopulationParameter> params = null;
	private List<String> symbols = new ArrayList<String>();
	
	/**
	 * Constructor
	 * @param cm_ Covariate Model
	 * @param lexer_ Converter Instance
	 */
	public CovariateBlockImpl(CovariateModel cm_, ILexer lexer_) {
		if (cm_ == null) throw new NullPointerException("The Covariate model cannot be null.");
		if (lexer_ == null) throw new NullPointerException("The lexer variable cannot be null.");
		
		cm = cm_;
		lexer = lexer_;
		
		params = cm.getListOfPopulationParameter();
		if (params != null) {
			for (PopulationParameter p : params) {
				if (p == null) continue;
				checkAssignment(p);
			}
		}
	}
	
	/**
	 * Associate a parameter with a covariate definition.
	 * @param covName Covariate Name
	 * @param parameterName Parameter Name
	 * @return boolean Success if named elements exist in the PharmML model.
	 */
	public boolean addParameterToCovariate(String covName, String parameterName) {
		if (covName == null || parameterName == null) return false;
		
		Accessor a = lexer.getAccessor();
		
		PharmMLRootType element = a.fetchElement(covName);
		if (!isCovariate(element)) return false;
		CovariateDefinition cov = (CovariateDefinition) element;
		
		element = a.fetchElement(parameterName);
		if (!isCommonParameter(element)) return false;
		CommonParameter p = (CommonParameter) element;
		
		CovariateParameterRef ref = null;
		if (!cov_param_map.containsKey(cov)) {
			ref = new CovariateParameterRef(cov);
			cov_param_map.put(cov, ref);
		} else
			ref = cov_param_map.get(cov);
		
		if (ref != null) return ref.addParameter(p);
		
		return false;
	}
	
	@Override
	public void buildTrees() {
		TreeMaker tm = lexer.getTreeMaker();
		
		if (params != null) {
			for (PopulationParameter p : params) {
				if (p == null) continue;
				BinaryTree bt = tm.newInstance(p);
				lexer.addStatement(p, bt);
			}
		}
		
		List<CovariateDefinition> covs = cm.getListOfCovariate();
		if (covs != null) {
			for (CovariateDefinition cov : covs) {
				if (cov == null) continue;
				
				lexer.addStatement(cov, tm.newInstance(cov));
				
				String symbId = cov.getSymbId();
				if (symbId == null) throw new NullPointerException("Covariate symbol ID is null");
				if (!symbols.contains(symbId)) symbols.add(symbId);
				
				if (cov.getCategorical() != null) {
					CategoricalCovariateRef cref = new CategoricalCovariateRef(lexer, cov, getName());
					categorical_covariates.add(cref);
					cref.buildTrees();
					lexer.addStatement(cref, tm.newInstance(cref));
					categorical_cov_map.put(cref.getName(), cref);
					categorical = true;
					categorical_cov_names.add(cov.getSymbId());
					
				} else if (cov.getContinuous() != null) {
					ContinuousCovariate continuous = cov.getContinuous();
					if (continuous.getInterpolation() != null) 
						throw new UnsupportedOperationException("covariate::interpolate not supported yet.");
					
					lexer.addStatement(continuous, tm.newInstance(continuous));
					lexer.updateNestedTrees();
										
					for (CovariateTransformation transform : continuous.getListOfTransformation()) 
						if (transform != null) lexer.addStatement(transform, tm.newInstance(transform));
					
					continuous_covariates.add(continuous);
				}
				covariates.add(cov);
			}
		}
	}
	
	/**
	 * Return a list of categorical covariate variable names
	 * @return java.util.List<String>
	 */
	public List<String> getCategoricalCovariateNames() { return categorical_cov_names; }
	
	/**
	 * Get a list of categorical covariate references.
	 * @return java.util.List<CategoricalCovariateRef>
	 */
	public List<CategoricalCovariateRef> getCategoricalCovariates() { return categorical_covariates; }
	
	@Override
	public List<String> getCategories(CovariateDefinition cov) { throw new UnsupportedOperationException(); }

	/**
	 * Get the continuous covariates from the model.
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modeldefn.ContinuousCovariate>
	 */
	public List<ContinuousCovariate> getContinuousCovariates() { return continuous_covariates; }
	
	/**
	 * Get a list of covariate to parameter references.
	 * @return Collection<CovariateParameterRef>
	 */
	public Collection<CovariateParameterRef> getCovariateParameterRefs() { return cov_param_map.values(); }
		
	/**
	 * Get a list of covariates in the model.
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition>
	 */
	public List<CovariateDefinition> getCovariates() { return covariates; }

	/**
	 * Get the source model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.CovariateModel
	 */
	public CovariateModel getModel() { return cm; }

	@Override
	public String getName() {
		String name = null;
		if (cm != null) name = cm.getBlkId();
		return name;
	}

	/**
	 * Get a list of parameters declared in the covariate block.
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter>
	 */
	public List<PopulationParameter> getParameters() {
		if (params == null) return new ArrayList<PopulationParameter>();
		else return params;
	}
	
	@Override
	public List<String> getSymbolIds() { return symbols; }
	
	
	@Override
	public boolean hasSymbolId(String name) {
		if (name != null) {
			if (symbols.contains(name)) return true;
		}
		
		return false;
	}
	
	/**
	 * Flag that the covariate model is a categorical.
	 * @return boolean
	 */
	public boolean isCategorical() { return categorical; }
	
	/**
	 * Flag that the 'named' covariate is categorical in nature.
	 * @return boolean
	 */
	public boolean isCategorical(String covariateName) {
		if (covariateName == null) return false;
		else return categorical_cov_map.containsKey(covariateName);
	}
}

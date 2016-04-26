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

import static crx.converter.engine.PharmMLTypeChecker.isIndividualParameter;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import crx.converter.engine.Accessor;
import crx.converter.engine.FixedParameter;
import crx.converter.engine.parts.SortableElement;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.steps.EstimationStep;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modellingsteps.Estimation;
import eu.ddmore.libpharmml.dom.modellingsteps.EstimationOpType;
import eu.ddmore.libpharmml.dom.modellingsteps.EstimationOperation;
import eu.ddmore.libpharmml.dom.modellingsteps.InitialEstimate;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;
import eu.ddmore.libpharmml.dom.modellingsteps.ToEstimate;

/**
 * A class the wraps an estimation code block.
 *
 */
public class EstimationStepImpl extends BaseStepImpl implements EstimationStep {
	private static double defaultParameterEstimateValue = 1.0;
	private static boolean useDefaultParameterEstimate = false;
	
	/**
	 * Check that the estimation operation is an individual parameter estimation.
	 * @param op
	 * @return boolean
	 * @see eu.ddmore.libpharmml.dom.modellingsteps.EstimationOpType#EST_INDIV
	 */
	public static boolean isIndividualParameterEstimate(EstimationOperation op) {
		if (op != null) {
			String type = op.getOpType();
			EstimationOpType value = EstimationOpType.fromValue(type);
			if (value != null) {
				return value.equals(EstimationOpType.EST_INDIV);
			}
		}
		
		return false;
	}
	
	/**
	 * Set the value for the default parameter estimate.<br/>
	 * Only assigned if not specified in a PharmML document.
	 * @param value Default parameter estimation value.
	 */
	public static void setDefaultParameterEstimateValue(double value) { defaultParameterEstimateValue = value; }
	
	/**
	 * Instruct a converter to use a default parameter estimate.<br/>
	 * Instruct the converter whether it can tolerate parameter calculation without an initial estimate.
	 * @param decision
	 */
	public static void setUseDefaultParameterEstimate(boolean decision) { useDefaultParameterEstimate = decision; }
	 
	private Estimation est = null;
	private HashMap<ParameterEstimate, Integer> estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private List<FixedParameter> fixed_parameters = new ArrayList<FixedParameter>();
	private HashMap<ParameterEstimate, Integer> indiv_estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private EstimationOperation [] operations = null;
	private List<ParameterEstimate> params_to_estimate = new ArrayList<ParameterEstimate>();
	
	/**
	 * Constructor
	 * @param est_ Estimation Step
	 * @param lexer_ Converter Instance
	 * @throws NullPointerException When constructor arguments NULL
	 * @throws IllegalStateException When estimation has no object identifier.
	 */
	public EstimationStepImpl(Estimation est_, ILexer lexer_) {
		if (est_ == null) throw new NullPointerException("The estimation step cannot be NULL");
		if (est_.getOid() == null)  throw new IllegalStateException("THe estimation step must have a OID");
		if (lexer_ == null) throw new NullPointerException("Lexer reference is NULL");
		
		est = est_;
		lexer = lexer_;
		a = lexer.getAccessor();
		exd = a.fetchExternalDataSet(est.getExternalDataSetReference());
	}
	
	private void buildOperationsArray() {
		// Look for operations associated with this step (if any)
		List<EstimationOperation> operations_ = est.getOperation();
		if (operations_ != null) {
			ArrayList<SortableElement> operation_refs = new ArrayList<SortableElement>();
			for (EstimationOperation operation : operations_) {
				if (operation == null) continue;
				if (operation.getOrder() == null) continue;
				operation_refs.add(new SortableElement(operation, operation.getOrder().intValue()));
			}
			Collections.sort(operation_refs);
			operations = new EstimationOperation[operation_refs.size()];
			for (int i = 0; i < operations.length; i++) 
				operations[i] = (EstimationOperation) operation_refs.get(i).getElement();
		} else 
			operations = new EstimationOperation[0];
	}
	
	private void buildParameterEstimateTrees() {
		Accessor a = lexer.getAccessor();
		
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) throw new NullPointerException("Model has no defined parameter block.");
		
		for (ParameterEstimate pe : params_to_estimate) {
			SymbolRef ref = pe.getSymbRef();
			if (ref == null) continue;
			
			PharmMLRootType element = a.fetchElement(ref);
			if (isPopulationParameter(element)) estimate_to_index.put(pe, pb.getParameterIndex(ref));
			else if (isIndividualParameter(element)) indiv_estimate_to_index.put(pe, pb.getParameterIndex(ref));
		}
	}
	
	@Override
	public void buildTrees() {
		if (est == null) throw new NullPointerException("Estimation step in model is NULL");
		
		categoriseParameterUsage();
		buildExternalDatasetMappings();
		buildOperationsArray();
		buildParameterEstimateTrees();
	}
	
	private void categoriseParameterUsage() {
		// Parameter estimates.
		TreeMaker tm = lexer.getTreeMaker();
		ToEstimate param_list_holder = est.getParametersToEstimate();
		if (param_list_holder != null) {
			List<ParameterEstimate> param_list = param_list_holder.getParameterEstimation();
			if (param_list != null) { 
				for (ParameterEstimate p : param_list) {
					if (p == null) continue;
					InitialEstimate ic = p.getInitialEstimate();
					if (ic ==  null && useDefaultParameterEstimate) ic = getDefaultInitialEstimate();

					boolean isFixed = false;
					if (ic != null) if (ic.isFixed() != null) isFixed = ic.isFixed();

					if (isFixed) {
						FixedParameter fp = new FixedParameter(p);
						fixed_parameters.add(fp);
						lexer.addStatement(fp, tm.newInstance(fp));
						lexer.updateNestedTrees(); 
					} else {					
						lexer.addStatement(p, tm.newInstance(p));
						lexer.updateNestedTrees(); 
						params_to_estimate.add(p);
					}
				}

				for (FixedParameter fp : fixed_parameters) param_list.remove(fp.pe);
			}
		}
	}
	
	private InitialEstimate getDefaultInitialEstimate() {
		InitialEstimate ic = new InitialEstimate();
		ic.setFixed(false);
		ic.setScalar(new RealValue(defaultParameterEstimateValue));
		
		return ic;
	}
	
	/**
	 * Retrieves the fixed parameter record associated with a parameter estimation.
	 * @param p Parameter Variable name
	 * @return FixedParameter
	 */
	public FixedParameter getFixedParameter(PopulationParameter p) {
		FixedParameter fp = null;
		
		if (p != null) {
			for (FixedParameter fp_ : fixed_parameters) {
				if (fp_ == null) continue;
				if (fp_.pe == null) continue;
				
				Object o = a.fetchElement(fp_.pe.getSymbRef());
				if (p.equals(o)) {
					fp = fp_;
					break;
				}
			}
		}
		
		return fp;
	}
	
	/**
	 * Get a list of fixed parameters.
	 * @return java.util.List<FixedParameter>
	 */
	public List<FixedParameter> getFixedParameters() {
		return fixed_parameters;
	}
	
	/**
	 * Get the index value for a number in a parameter estimate vector.
	 * @param pe Parameter Estimate
	 * @return java.lang.Integer
	 */
	public Integer getIndividualParameterIndex(ParameterEstimate pe) {
		Integer idx = -1;
		
		if (pe != null) {
			if (indiv_estimate_to_index.containsKey(pe)) idx = indiv_estimate_to_index.get(pe);
		}
		
		return idx;
	}
		
	@Override
	public String getName() {
		String blkId = "__RUBBISH_DEFAULT_VALUE_";
		if (est != null) {
			if (est.getOid()  != null) blkId = est.getOid();
		}
		return blkId;
	}
	
	/**
	 * List of operations associated with this estimation step.
	 * @return eu.ddmore.libpharmml.dom.modellingsteps.EstimationOperation[]
	 */
	public EstimationOperation [] getOperations() { return operations; }
	
	/**
	 * Get the parameter estimate objedt associated with a parameter object.
	 * @param p Parameter
	 * @return eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate
	 */
	public ParameterEstimate getParameterEstimate(PopulationParameter p) {
		ParameterEstimate pe = null;
		
		if (p != null) {
			for (ParameterEstimate pe_ : getParametersToEstimate()) {
				if (pe_ == null) continue;
				PharmMLRootType element = a.fetchElement(pe_.getSymbRef());
				if (element == null) continue;
				if (element.equals(p)) {
					pe = pe_;
					break;
				}
			}
		}
		
		return pe;
	}
	
	/**
	 * Get the index number of a parameter in an estimation vector.
	 * @param pe Parameter Estimate
	 * @return java.lang.Integer
	 */
	public Integer getParameterIndex(ParameterEstimate pe) {
		Integer idx = -1;
		
		if (pe != null) {
			if (estimate_to_index.containsKey(pe)) idx = estimate_to_index.get(pe);
		}
		
		return idx;
	}
	
	/**
	 * Get the list of parameter estimates
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate>
	 */
	public List<ParameterEstimate> getParametersToEstimate() { return params_to_estimate; }

	/**
	 * Get the estimation step.
	 * @return eu.ddmore.libpharmml.dom.modellingsteps.Estimation
	 */
	public Estimation getStep() { return est; }
	
	@Override
	public List<String> getSymbolIds() { return new ArrayList<String>(); }

	/**
	 * Flag if the estimation has fixed parameters.
	 * @return boolean
	 */
	public boolean hasFixedParameters() { return fixed_parameters.size() > 0; }
	
	/**
	 * If the estimation has parameters to estimate.
	 * @return boolean
	 */
	public boolean hasParametersToEstimate() { return params_to_estimate.size() > 0; }
	
	/**
	 * Check if the estimation has simple parameters to estimate.
	 * @return boolean
	 */
	public boolean hasSimpleParametersToEstimate() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public boolean hasSymbolId(String name) {
		return false;
	}
	
	/**
	 * Flag if the estimation is constrained.
	 * @return boolean
	 */
	public boolean isConstrained() { throw new UnsupportedOperationException(); }
	
	/**
	 * Check is a given parameter is a fixed parameter in an estimation.
	 * @param p
	 * @return boolean
	 */
	public boolean isFixedParameter(PopulationParameter p) {
		boolean isFixed = false;
		if (p != null) {
			for (FixedParameter fixed_parameter : fixed_parameters) {
				if (fixed_parameter == null) continue;
				if (fixed_parameter.pe != null) {
					PharmMLRootType element = a.fetchElement(fixed_parameter.pe.getSymbRef());
					if (element != null) {
						if (p.equals(element)) {
							isFixed = true;
							break;
						}
					}
				}
			}
		}
		
		return isFixed;
	}
	
	/**
	 * Update the parameter estimation order to that expressed in the parameter model.
	 */
	public void update() {
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		List<ParameterEstimate> pes = new ArrayList<ParameterEstimate>();
		pes.addAll(params_to_estimate);
		params_to_estimate.clear();
		estimate_to_index.clear();
		
		Accessor a = lexer.getAccessor();
		int idx = 0;
		for (PopulationParameter p : pb.getParameters()) {
			if (p == null) continue;
			ParameterEstimate pe_found = null;
			for (ParameterEstimate pe : pes) {
				if (pe == null) continue;
				SymbolRef ref = pe.getSymbRef();
				if (ref != null) {
					Object o = a.fetchElement(ref);
					if (p.equals(o)) {
						pe_found = pe;
						break;
					}
				}
			}
			
			if (pe_found != null) {
				params_to_estimate.add(pe_found);
				estimate_to_index.put(pe_found, idx++);
				pes.remove(pe_found);
			}
		}
		
		for (ParameterEstimate pe : pes) {
			if (pe == null) continue;
			params_to_estimate.add(pe);
			estimate_to_index.put(pe, idx++);
		}
	}
	
	@Override
	public String getToolName() { throw new UnsupportedOperationException(); }
}

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
import static eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOpType.EVALUATION;
import static eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOpType.OPTIMISATION;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.engine.FixedParameter;
import crx.converter.engine.parts.SortableElement;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.steps.OptimalDesignStep_;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modellingsteps.InitialEstimate;
import eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOpType;
import eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOperation;
import eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignStep;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;
import eu.ddmore.libpharmml.dom.modellingsteps.ToEstimate;

/**
 * PFIM specific implementation of an optimal design step wrapper.
 */
public class OptimalDesignStepImpl extends BaseStepImpl implements OptimalDesignStep_ {
	private static double defaultParameterEstimateValue = 1.0;
	private static boolean useDefaultParameterEstimate = false;
	
	/**
	 * Set the value for the default parameter estimate.<br/>
	 * Only assigned if not specified in a PharmML document.
	 * @param value Default parameter estimation value.
	 */
	public static void setDefaultParameterEstimateValue(double value) { defaultParameterEstimateValue = value; }
	
	private Map<ParameterEstimate, Integer> estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private boolean evaluation = false;
	private List<FixedParameter> fixed_parameters = new ArrayList<FixedParameter>();
	private Map<ParameterEstimate, Integer> indiv_estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private OptimalDesignOperation [] operations = null;
	private boolean optimisation = false;
	private List<ParameterEstimate> params_to_estimate = new ArrayList<ParameterEstimate>();
	private OptimalDesignStep step = null;
	
	/**
	 * Constructor
	 * @param step Optimal Design Step
	 * @param lexer Converter/Lexer handle
	 */
	public OptimalDesignStepImpl(OptimalDesignStep step_, ILexer lexer_) {
		if (step_ == null) throw new NullPointerException("The optimal design step cannot be NULL");
		if (lexer_ == null) throw new NullPointerException("Lexer/Converter cannot be NULL");
		
		step = step_;
		lexer = lexer_;
		a = lexer_.getAccessor();
	}
	
	private void buildOperationsArray() {
		// Look for operations associated with this step (if any)
		List<OptimalDesignOperation> operations_ = step.getListOfOperation();
		if (operations_ != null) {
			ArrayList<SortableElement> operation_refs = new ArrayList<SortableElement>();
			for (OptimalDesignOperation operation : operations_) {
				if (operation == null) continue;
				if (operation.getOrder() == null) continue;
				operation_refs.add(new SortableElement(operation, operation.getOrder().intValue()));
			}
			Collections.sort(operation_refs);
			operations = new OptimalDesignOperation[operation_refs.size()];
			for (int i = 0; i < operations.length; i++) 
				operations[i] = (OptimalDesignOperation) operation_refs.get(i).getElement();
		} else 
			operations = new OptimalDesignOperation[0];
	}
	
	private void buildParameterEstimateTrees() {
		Accessor a = lexer.getAccessor();
		TreeMaker tm = lexer.getTreeMaker();
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) throw new NullPointerException("Model has no defined parameter block.");
		
		for (ParameterEstimate pe : params_to_estimate) {
			if (pe == null) continue;
			SymbolRef ref = pe.getSymbRef();
			if (ref == null) continue;
			
			PharmMLRootType element = a.fetchElement(ref);
			if (isPopulationParameter(element)) {
				estimate_to_index.put(pe, pb.getParameterIndex(ref));
			} else if (isIndividualParameter(element)) indiv_estimate_to_index.put(pe, pb.getParameterIndex(ref));
			
			lexer.addStatement(pe, tm.newInstance(pe));
		}
	}
	
	@Override
	public void buildTrees() {
		categoriseParameterUsage();
		buildOperationsArray();
		setTaskType();
		buildParameterEstimateTrees();
	}
	
	private void categoriseParameterUsage() {
		// Parameter estimates.
		TreeMaker tm = lexer.getTreeMaker();
		ToEstimate param_list_holder = step.getParametersToEstimate();
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
					}
					
					params_to_estimate.add(p);
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
	public List<FixedParameter> getFixedParameters() { return fixed_parameters; }

	@Override
	public String getName() { return step.getOid(); }

	/**
	 * Get the operations array bound to the OD step.
	 * @return
	 */
	public OptimalDesignOperation[] getOperations() { return operations; }
	
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
	 * Get the source PharmML step for the optimal design.
	 * @return OptimalDesignStep
	 */
	public OptimalDesignStep getStep() { return step; }

	@Override
	public List<String> getSymbolIds() { return new ArrayList<String>(); }

	@Override
	public String getToolName() { return ""; }

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

	@Override
	public boolean hasSymbolId(String name) { return false; }

	/**
	 * Flag if the parameter estimation is constrained.
	 * @return boolean
	 */
	public boolean isConstrained() {
		for (ParameterEstimate pe : params_to_estimate) {
			if (pe == null) continue;
			if (pe.getLowerBound() != null || pe.getUpperBound() != null) return true;
		}
		
		return false;
	}

	/**
	 * Flag if an evaluation task.
	 * @return boolean
	 */
	public boolean isEvaluation() { return evaluation; }

	@Override
	public boolean isFixed(IndividualParameter ip) {
		if (ip == null) return false;
		
		if (ip.getStructuredModel() != null) return false;
		else if (ip.getAssign() != null) return true;
		
		return false;
	}
	
	/**
	 * Flag if an optimisation task.
	 * @return boolean
	 */
	public boolean isOptimisation() { return optimisation; }
	
	private void setTaskType() {
		if (operations == null) return;
		
		for (OptimalDesignOperation op : operations) {
			if (op == null) continue;
			OptimalDesignOpType type = op.getOpType();
			if (type == null) continue;
			if (EVALUATION.equals(type)) evaluation = true;
			else if (OPTIMISATION.equals(type)) optimisation = true;
		}
	}
	
	@Override
	public String toString() { return step.getOid(); }
}

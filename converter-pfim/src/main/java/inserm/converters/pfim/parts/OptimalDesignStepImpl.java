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
import static crx.converter.engine.PharmMLTypeChecker.isRandomVariable;
import static eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOpType.EVALUATION;
import static eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignOpType.OPTIMISATION;
import static inserm.converters.pfim.SettingLabel.ALGORITHM;
import inserm.converters.pfim.OptimisationAlgorithm;
import inserm.converters.pfim.Parser;
import inserm.converters.pfim.SettingReader;

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
import crx.converter.spi.blocks.VariabilityBlock;
import crx.converter.spi.steps.OptimalDesignStep_;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.LevelReference;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomEffect;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.ParentLevel;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel;
import eu.ddmore.libpharmml.dom.modeldefn.VariabilityLevelDefinition;
import eu.ddmore.libpharmml.dom.modellingsteps.Algorithm;
import eu.ddmore.libpharmml.dom.modellingsteps.InitialEstimate;
import eu.ddmore.libpharmml.dom.modellingsteps.OperationProperty;
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
	
	private String algorithm = null;
	private Map<ParameterEstimate, Integer> estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private boolean evaluation = false;
	private List<FixedParameter> fixed_parameters = new ArrayList<FixedParameter>();
	private Map<IndividualParameter, List<ParameterRandomVariable>> gamma_map = new HashMap<IndividualParameter, List<ParameterRandomVariable>>();
	private Map<ParameterEstimate, Integer> indiv_estimate_to_index = new HashMap<ParameterEstimate, Integer>();
	private Map<IndividualParameter, List<ParameterRandomVariable>> omega_map = new HashMap<IndividualParameter, List<ParameterRandomVariable>>();
	private OptimalDesignOperation [] operations = null;
	private boolean optimisation = false;
	private List<ParameterEstimate> params_to_estimate = new ArrayList<ParameterEstimate>();
	private SettingReader sr = null;
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
	
	private void buildGammaMap() {
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		VariabilityBlock vb = getParameterVariabilityScope();
		if (vb == null) return;
		List<VariabilityLevelDefinition> levels = vb.getModel().getLevel();
		if (levels == null) return;
		if (levels.size() > 2) throw new IllegalStateException("More than the expected variability levels");
		
		String id_level = null, occassion_level = null;
		for (VariabilityLevelDefinition level : levels) {
			if (level == null) continue;
			if (level.getParentLevel() == null) {
				id_level = level.getSymbId();
				break;
			}
		}
		if (id_level == null) return;
		for (VariabilityLevelDefinition level : levels) {
			if (level == null) continue;
			ParentLevel parent = level.getParentLevel();
			if (parent != null) {
				SymbolRef ref = parent.getSymbRef();
				if (ref != null) {
					if (id_level.equalsIgnoreCase(ref.getSymbIdRef())) {
						occassion_level = level.getSymbId();
						break;
					}
				}
			}
		}
		if (occassion_level == null) return;
		
		Accessor a = lexer.getAccessor();
		List<IndividualParameter> ips = pb.getIndividualParameters();
		for (IndividualParameter ip : ips) {
			if (ip == null) throw new NullPointerException("An individual parameter is NULL");
			if (!isFixed(ip)) {
				List<ParameterRandomVariable> gammas = new ArrayList<ParameterRandomVariable>();
				StructuredModel soe = ip.getStructuredModel();
				if (soe == null) continue;
				List<ParameterRandomEffect> rves = soe.getListOfRandomEffects();
				if (rves == null) continue;
				if (rves.size() > 2) throw new IllegalStateException("More than expected random effects in individual variable (symbId='" + ip.getSymbId() + "')");
				for (ParameterRandomEffect rve : rves) {
					if (rve == null) continue;
					if (rve.getSymbRef().isEmpty()) continue;
					SymbolRef ref = rve.getSymbRef().get(0);
					if (ref != null) {
						PharmMLRootType element = a.fetchElement(ref);
						if (isRandomVariable(element)) {
							ParameterRandomVariable rv = (ParameterRandomVariable) element;
							for (LevelReference lref : rv.getListOfVariabilityReference()) {
								if (lref == null) continue;
								ref = lref.getSymbRef();
								if (ref != null) {
									if (occassion_level.equalsIgnoreCase(ref.getSymbIdRef())) gammas.add(rv);
								}
							}
						}
					}
				}
				
				if (gammas.size() > 0) gamma_map.put(ip, gammas);
			}
		}
	}
	
	// PFIM only support 1 omega per individual variable.
	private void buildOmegaMap() {
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		VariabilityBlock vb = getParameterVariabilityScope();
		if (vb == null) return;
		List<VariabilityLevelDefinition> levels = vb.getModel().getLevel();
		if (levels == null) return;
		if (levels.size() > 2) throw new IllegalStateException("More than the expected variability levels");
		
		String id_level = null;
		for (VariabilityLevelDefinition level : levels) {
			if (level == null) continue;
			if (level.getParentLevel() == null) {
				id_level = level.getSymbId();
				break;
			}
		}
		if (id_level == null) return;
		
		Accessor a = lexer.getAccessor();
		List<IndividualParameter> ips = pb.getIndividualParameters();
		for (IndividualParameter ip : ips) {
			if (ip == null) throw new NullPointerException("An individual parameter is NULL");
			if (!isFixed(ip)) {
				List<ParameterRandomVariable> omegas = new ArrayList<ParameterRandomVariable>();
				StructuredModel soe = ip.getStructuredModel();
				if (soe == null) continue;
				List<ParameterRandomEffect> rves = soe.getListOfRandomEffects();
				if (rves == null) continue;
				if (rves.size() > 2) throw new IllegalStateException("More than expected random effects in individual variable (symbId='" + ip.getSymbId() + "')");
				for (ParameterRandomEffect rve : rves) {
					if (rve == null) continue;
					if (rve.getSymbRef().isEmpty()) continue;
					SymbolRef ref = rve.getSymbRef().get(0);
					if (ref != null) {
						PharmMLRootType element = a.fetchElement(ref);
						if (isRandomVariable(element)) {
							ParameterRandomVariable rv = (ParameterRandomVariable) element;
							for (LevelReference lref : rv.getListOfVariabilityReference()) {
								if (lref == null) continue;
								ref = lref.getSymbRef();
								if (ref != null) {
									if (id_level.equalsIgnoreCase(ref.getSymbIdRef())) omegas.add(rv);
								}
							}
						}
					}
				}
				
				if (omegas.size() > 0) omega_map.put(ip, omegas);
			}
		}
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
		buildOmegaMap();
		buildGammaMap();
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

	@Override
	public String getAlgorithm() { return algorithm; }

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
	public List<ParameterRandomVariable> getGammas(IndividualParameter ip) {
		if (ip != null) {
			if (gamma_map.containsKey(ip)) return gamma_map.get(ip); 
		}
		
		return null;
	}

	@Override
	public String getName() { return step.getOid(); }

	@Override
	public List<ParameterRandomVariable> getOmegas(IndividualParameter ip) {
		if (ip != null) {
			if (omega_map.containsKey(ip)) return omega_map.get(ip); 
		}
		
		return null;
	}

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

	private VariabilityBlock getParameterVariabilityScope() {
		VariabilityBlock vb = null;
		
		List<VariabilityBlock> vbs = lexer.getScriptDefinition().getVariabilityBlocks();
		for (VariabilityBlock vb_ : vbs) {
			if (vb_ == null) continue;
			if (vb_.isParameterVariability()) {
				vb = vb_;
				break;
			}
		}
		
		return vb;
	}

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

	// Read the algorithm and other linked optimal design settings linked to PFIM.
	private void readAlogrithm(OptimalDesignOperation op) {
		if (op == null) return;
		
		Algorithm algo = op.getAlgorithm();
		if (algo == null) return;
		List<OperationProperty> properties =  algo.getProperty();
		if (properties == null) return;
		if (properties.isEmpty()) return;
		
		sr = new SettingReader();
		sr.setLexer(lexer);
		sr.setParser(lexer.getParser());
		sr.setProperties(properties);
		sr.readSettings();
		
		Parser p = (Parser) lexer.getParser();
		p.register(sr);
		
		if(!sr.hasValue(ALGORITHM.toString())) return;
		String definition = sr.getValue(ALGORITHM.toString());
		if (definition == null) throw new NullPointerException("OD algorithm definition not specified.");
		if (definition.isEmpty()) throw new NullPointerException("OD algorithm definition is a zero-length string.");
		
		definition = definition.toUpperCase();
		if (!OptimisationAlgorithm.contains(definition)) 
			throw new UnsupportedOperationException("The specified algorithm is not supported by PFIM (algo='" + definition + "')");
		else
			algorithm = definition;
	}

	private void setTaskType() {
		if (operations == null) return;
		
		for (OptimalDesignOperation op : operations) {
			if (op == null) continue;
			OptimalDesignOpType type = op.getOpType();
			if (type == null) continue;
			if (EVALUATION.equals(type)) evaluation = true;
			else if (OPTIMISATION.equals(type)) optimisation = true;
			
			readAlogrithm(op);
		}
	}

	@Override
	public String toString() { return step.getOid(); }
}

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
import static crx.converter.engine.PharmMLTypeChecker.isCorrelation;
import static crx.converter.engine.PharmMLTypeChecker.isCovariate;
import static crx.converter.engine.PharmMLTypeChecker.isDerivative;
import static crx.converter.engine.PharmMLTypeChecker.isFunctionCall;
import static crx.converter.engine.PharmMLTypeChecker.isGeneralError;
import static crx.converter.engine.PharmMLTypeChecker.isIndividualParameter;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isObservationError;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;
import static crx.converter.engine.PharmMLTypeChecker.isRandomVariable;
import static crx.converter.engine.PharmMLTypeChecker.isStructuredError;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;

import java.util.ArrayList;
import java.util.List;

import crx.converter.engine.Accessor;
import crx.converter.engine.common.CorrelationRef;
import crx.converter.engine.common.ObservationParameter;
import crx.converter.engine.common.SimulationOutput;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.tree.NestedTreeRef;
import crx.converter.tree.Node;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLElement;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.StandardAssignable;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType.FunctionArgument;
import eu.ddmore.libpharmml.dom.modeldefn.CategoricalData;
import eu.ddmore.libpharmml.dom.modeldefn.CategoricalPMF;
import eu.ddmore.libpharmml.dom.modeldefn.Censoring;
import eu.ddmore.libpharmml.dom.modeldefn.CommonDiscreteVariable;
import eu.ddmore.libpharmml.dom.modeldefn.CommonObservationModel;
import eu.ddmore.libpharmml.dom.modeldefn.CommonParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ContinuousObservationModel;
import eu.ddmore.libpharmml.dom.modeldefn.Correlation;
import eu.ddmore.libpharmml.dom.modeldefn.CountData;
import eu.ddmore.libpharmml.dom.modeldefn.CountPMF;
import eu.ddmore.libpharmml.dom.modeldefn.Dependance;
import eu.ddmore.libpharmml.dom.modeldefn.Discrete;
import eu.ddmore.libpharmml.dom.modeldefn.DiscreteDataParameter;
import eu.ddmore.libpharmml.dom.modeldefn.GeneralObsError;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ObservationError;
import eu.ddmore.libpharmml.dom.modeldefn.ObservationModel;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ProbabilityAssignment;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError;
import eu.ddmore.libpharmml.dom.modeldefn.TTEFunction;
import eu.ddmore.libpharmml.dom.modeldefn.TimeToEventData;

/**
 * The observation block of a PharmML model.<br/>
 *
 */
public class ObservationBlockImpl extends BaseRandomVariableBlockImpl implements ObservationBlock {
	private Accessor a = null;
	private CategoricalData categorical_data = null;
	private CountData count_data = null;
	private Dependance dependence = null;
	private Discrete discrete = null;
	private ObservationError error_model = null;
	private List<NestedTreeRef> error_trees = new ArrayList<NestedTreeRef>();
	private ArrayList<String> func_names = new ArrayList<String>();
	private ArrayList<String> model_block_ids = new ArrayList<String>();
	private ObservationModel ob = null;
	private ArrayList<ObservationParameter> param_refs = new ArrayList<ObservationParameter>();
	private ArrayList<CommonParameter> params = new ArrayList<CommonParameter>();
	private ArrayList<ParameterRandomVariable> random_variables = new ArrayList<ParameterRandomVariable>();
	private ArrayList<SimulationOutput> simulation_outputs = new ArrayList<SimulationOutput>();
	private TreeMaker tm = null;
	private TimeToEventData tte_data = null; 
	
	/**
	 * Constructor
	 * @param ob_ Observation Model
	 * @param lexer_ Converter Reference
	 */
	public ObservationBlockImpl(ObservationModel ob_, ILexer lexer_) {
		if (lexer_ == null) throw new NullPointerException("The lexer cannot be NULL.");
		if (ob_ == null) throw new NullPointerException("The observation model type cannot be NULL.");
		
		ob = ob_;
		lexer = lexer_;
		a = lexer.getAccessor();
		tm = lexer.getTreeMaker();

		initErrorModel();
	}
	
	private void buildCategoricalData() {
		if (categorical_data == null) return; 
		
		CommonDiscreteVariable cv = categorical_data.getCategoryVariable();
		if (cv != null) {
			lexer.addStatement(cv, tm.newInstance(cv));
			lexer.updateNestedTrees();
		}
		
		dependence = categorical_data.getDependance();
		
		List<CommonDiscreteVariable> categories = categorical_data.getListOfCategories();
		if (categories != null) {
			for (CommonDiscreteVariable category : categories) {
				if (category == null) continue;
				lexer.addStatement(category, tm.newInstance(category));
				lexer.updateNestedTrees();
			}
		}
		
		List<CategoricalPMF> pmfs = categorical_data.getListOfPMF();
		if (pmfs != null) {
			for (CategoricalPMF pmf : pmfs) {
				if (pmf == null) continue;
				lexer.addStatement(pmf, tm.newInstance(pmf));
				lexer.updateNestedTrees();
			}
		}
		
		List<CommonDiscreteVariable> previousStateVariables = categorical_data.getListOfPreviousStateVariable();
		if (previousStateVariables != null) {
			for (CommonDiscreteVariable previousStateVariable : previousStateVariables) {
				if (previousStateVariable == null) continue;
				lexer.addStatement(previousStateVariable, tm.newInstance(previousStateVariable));
				lexer.updateNestedTrees();
			}
		}
		
		List<ProbabilityAssignment> pas = categorical_data.getListOfProbabilityAssignment();
		if (pas != null) {
			for (ProbabilityAssignment pa : pas) {
				if (pa == null) continue;
				lexer.addStatement(pa, tm.newInstance(pa));
				lexer.updateNestedTrees();
			}
		}
		
		buildDiscreteCommon(categorical_data);
	}
	
	private void buildCorrelations() { for (CorrelationRef corr : correlations) buildTrees(corr); }
	
	private void buildCountData() {
		if (count_data == null) return;
		
		CommonDiscreteVariable count_variable = count_data.getCountVariable();
		if (count_variable != null) {
			lexer.addStatement(count_variable, tm.newInstance(count_variable));
			lexer.updateNestedTrees();
		}
		
		dependence = count_data.getDependance();
		
		DiscreteDataParameter dp = count_data.getDispersionParameter();
		if (dp != null) {
			lexer.addStatement(dp, tm.newInstance(dp));
			lexer.updateNestedTrees();
		}
		
		List<DiscreteDataParameter> ips = count_data.getListOfIntensityParameter();
		if (ips != null) {
			for (DiscreteDataParameter ip : ips) {
				if (ip == null) continue;
				lexer.addStatement(ip, tm.newInstance(ip));
				lexer.updateNestedTrees();
			}
		}
		
		List<CountPMF> pmfs = count_data.getListOfPMF();
		if (pmfs != null) {
			for (CountPMF pmf : pmfs) {
				if (pmf == null) continue;
				lexer.addStatement(pmf, tm.newInstance(pmf));
				lexer.updateNestedTrees();
			}
		}
		
		List<CommonDiscreteVariable> previous_counts = count_data.getListOfPreviousCountVariable();
		if (previous_counts != null) {
			for (CommonDiscreteVariable previous_count : previous_counts) {
				if (previous_count == null) continue;
				lexer.addStatement(previous_count, tm.newInstance(previous_count));
				lexer.updateNestedTrees();
			}
		}
		
		DiscreteDataParameter p = count_data.getMixtureProbabilityParameter();
		if (p != null) {
			lexer.addStatement(p, tm.newInstance(p));
			lexer.updateNestedTrees();
		}
		
		p = count_data.getOverDispersionParameter();
		if (p != null) {
			lexer.addStatement(p, tm.newInstance(p));
			lexer.updateNestedTrees();
		}
		
		p = count_data.getZeroProbabilityParameter();
		if (p != null) {
			lexer.addStatement(p, tm.newInstance(p));
			lexer.updateNestedTrees();
		}
		
		buildDiscreteCommon(count_data);
	}
	
	private void buildDiscreteCommon(CommonObservationModel com) {
		if (com == null) return;
		
		List<VariableDefinition> vs = new ArrayList<VariableDefinition>();
		for (PharmMLElement o : com.getListOfObservationModelElement()) {
			if (isCommonParameter(o)) params.add((CommonParameter) o);
			else if (isLocalVariable(o)) vs.add((VariableDefinition) o);
			else if (isCorrelation(o)) register((Correlation) o);
		}
	}
	
	private void buildErrorModel() {
		if (isStructuredError(error_model)) { 
			StructuredObsError goe = (StructuredObsError) error_model;
			error_trees.add(new NestedTreeRef(goe, tm.newInstance(goe)));
			
			StructuredObsError.ErrorModel error = goe.getErrorModel();
			if (error == null) throw new NullPointerException("Gaussian erorr model not specified.");
			error_trees.add(new NestedTreeRef(error, tm.newInstance(error.getAssign())));
			
			StructuredObsError.Output output = goe.getOutput();
			if (output == null) throw new NullPointerException("Gaussian erorr model, output variable not specified.");
			if (output.getSymbRef() == null) throw new NullPointerException("Gaussian Output variable not specified (symbId='" +  goe.getSymbId() +"')");
				
			StructuredObsError.ResidualError residual_error = goe.getResidualError();
			if (residual_error == null) throw new NullPointerException("Gaussian erorr model, residual variable not specified.");	
			if (residual_error.getSymbRef() == null) throw new NullPointerException("Gaussian Residual error variable not specified (symbId='" +  goe.getSymbId() +"')");
			
			// Add the referenced output variable to the simulation output list.
			PharmMLRootType element = a.fetchElement(output.getSymbRef());
			if (element == null) throw new NullPointerException("Element referenced by error model is NULL.");
			createSimulationOutput(element);
		} else if (isGeneralError(error_model)) {
			GeneralObsError goe = (GeneralObsError) error_model;
			error_trees.add(new NestedTreeRef(goe, tm.newInstance(goe)));
			error_trees.add(new NestedTreeRef(goe.getAssign(), tm.newInstance(goe.getAssign())));
		}
	}
	
	private void buildNestedTrees() {
		ArrayList<String> func_names = new ArrayList<String>();
		for (NestedTreeRef nref : error_trees) {
			for (Node node : nref.bt.nodes) {
				if (node == null) continue;
				
				if (isSymbolReference(node.data)) { 
					PharmMLRootType element = a.fetchElement((SymbolRef) node.data);
					if (element == null) continue;
					if (isDerivative(element) || isLocalVariable(element)) {
						SimulationOutput simulation_output = createSimulationOutput(element);
						node.data = simulation_output;
					}
					else if (isRandomVariable(element) || isPopulationParameter(element) || isIndividualParameter(element) || isCovariate(element)) continue;
					else throw new UnsupportedOperationException("Class not supported in 'observation' error model (" + element + ")");	
				}  
				else if (isFunctionCall(node.data)) {
					FunctionCallType call = (FunctionCallType) node.data;
					SymbolRef ref = call.getSymbRef();
					if (ref != null) {
						String func_name = ref.getSymbIdRef();
						if (func_name != null) {
							if (!func_names.contains(func_name)) func_names.add(func_name);
						}
					}
					
					List<FunctionArgument> args = call.getListOfFunctionArgument();
					if (args != null) {
						if (!args.isEmpty()) {
							for (FunctionArgument arg : args) {
								if (arg == null) continue;
								SymbolRef arg_ref =  arg.getSymbRef();
								if (arg_ref == null) continue;
								PharmMLRootType element = a.fetchElement(arg_ref);
								if (isDerivative(element) || isLocalVariable(element)) {
									SimulationOutput simulation_output = createSimulationOutput(element);
									node.data = simulation_output;
								}
							}
						}
					}
				}	
			}
			
			lexer.addStatement(nref);
		}
	}
	
	private void buildParameters() {
		for (CommonParameter param : params) {
			if (param == null) continue;
			if (isPopulationParameter(param)) {
				PopulationParameter sp = (PopulationParameter) param;
				ObservationParameter param_ref = new ObservationParameter(this, sp);
				param_refs.add(param_ref);
				lexer.addStatement(param_ref, tm.newInstance(sp));
			} else if (isRandomVariable(param)) {
				ParameterRandomVariable rv = (ParameterRandomVariable) param;
				lexer.addStatement(rv, tm.newInstance(rv));
				lexer.updateNestedTrees();
			}
		}
	}
	
	private void buildTimeToEventData() {
		if (tte_data == null) return;
		
		CommonDiscreteVariable ev = tte_data.getEventVariable();
		if (ev != null) {
			lexer.addStatement(ev, tm.newInstance(ev));
			lexer.updateNestedTrees();
		}
		
		List<Censoring> censorings = tte_data.getListOfCensoring();
		if (censorings != null) {
			for (Censoring censoring : censorings) {
				if (censoring == null) continue;
				lexer.addStatement(censoring, tm.newInstance(censoring));
				lexer.updateNestedTrees();
			}
		}
		
		List<StandardAssignable> features = tte_data.getListOfMaximumNumberEvents();
		if (features != null) {
			for (StandardAssignable feature : features) {
				if (feature == null) continue;
				lexer.addStatement(feature, tm.newInstance(feature));
				lexer.updateNestedTrees();
			}
		}
		
		List<TTEFunction> hazard_funcs = tte_data.getListOfHazardFunction();
		if (hazard_funcs != null) {
			for (TTEFunction hazard_func : hazard_funcs) {
				if (hazard_func == null) continue;
				lexer.addStatement(hazard_func, tm.newInstance(hazard_func));
				lexer.updateNestedTrees();
			}
		}
		
		List<TTEFunction> survival_funcs = tte_data.getListOfSurvivalFunction();
		if (survival_funcs != null) {
			for (TTEFunction survival_func : survival_funcs) {
				if (survival_funcs == null) continue;
				lexer.addStatement(survival_func, tm.newInstance(survival_func));
				lexer.updateNestedTrees();
			}
		}
		
		buildDiscreteCommon(count_data);		
	}
	
	@Override
	public void buildTrees() {
		buildCountData();
		buildCategoricalData();
		buildTimeToEventData();
		buildParameters();
		buildCorrelations();
		buildIndividualParameterTrees();
		buildErrorModel();
		buildNestedTrees();
	}
	
	/**
	 * Check if the error model contains a specific model element.
	 * @param v Model Element
	 * @return boolean
	 */
	public boolean contains(PharmMLRootType v) {
		if (isPopulationParameter(v)) {
			if (params.contains(v)) return true;
		} else {
			if (isObservationError(v)) {
				if (error_model != null) {
					if (error_model.equals(v)) return true;
				}
			}
		}
			
		return false;
	}
	
	private SimulationOutput createSimulationOutput(PharmMLRootType element) {
		if (element == null) return null;
		SimulationOutput output = null;
		for (SimulationOutput simulation_output : simulation_outputs) {
			if (simulation_output == null) continue;
			if (element.equals(simulation_output.v))  {
				output = simulation_output;
				break;
			}
		}
		
		if (output == null) {
			String blkId = a.getBlockId(element);
			if (blkId != null) {
				if (!model_block_ids.contains(blkId)) model_block_ids.add(blkId);
			}
			
			output = new SimulationOutput(element);
			simulation_outputs.add(output);
		}
		
		return output;
	}
	
	/**
	 * Get the categorical data structure of a discrete model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.CategoricalData
	 */
	public CategoricalData getCategoricalData() { return categorical_data; }
	
	/**
	 * Get the source count data structure associated with a discrete model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.CountData
	 */
	public CountData getCountData() { return count_data; }
		
	/**
	 * Get the dependance bound to a discrete model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.Dependance
	 */
	public Dependance getDependance() { return dependence; }
	
	/**
	 * Get a list of function names used in the error model.<br/>
	 * Cross-reference terms against a function library.
	 * @return boolean
	 */
	public List<String> getErrorFunctionNames() { return func_names; }
	
	/**
	 * Get the source observation model for an enclosing observation block.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.ObservationModel
	 */
	public ObservationModel getModel() { return ob; }
	
	@Override
	public String getName() {
		String blkId = "__RUBBISH_DEFAULT_VALUE_";
		if (ob != null) {
			if (ob.getBlkId() != null) blkId = ob.getBlkId();
		}
		return blkId;
	}
	
	@Override
	public ObservationError getObservationError() { return error_model; }
	
	@Override
	public ObservationParameter getObservationParameter(PopulationParameter p) {
		ObservationParameter op = null;
		
		if (p != null) {
			for (ObservationParameter param_ref : param_refs) {
				if (param_ref == null) continue;
				if (param_ref.param != null) {		
					if (param_ref.param.equals(p)) {
						op = param_ref;
						break;
					}
				}
			}
		}
		
		return op;
	}
	
	/**
	 * List of observation-model scoped parameters.
	 * @return java.util.List<ObservationParameter>
	 */
	public List<ObservationParameter> getObservationParameters() {
		return param_refs;
	}
	
	/**
	 * Get list of random variables declared in the error model.
	 * @return java.util.List<ParameterRandomVariable>
	 */
	public List<ParameterRandomVariable> getRandomVariables() {
		return random_variables;
	}
	
	/**
	 * Get the list of simulation outputs referenced by the error model.
	 * @return java.util.List<SimulationOutput>
	 */
	public List<SimulationOutput> getSimulationOutputs() {
		return simulation_outputs;
	}
	
	@Override
	public List<String> getSymbolIds() {
		ArrayList<String> ids = new ArrayList<String>();
		
		for (CommonParameter p : params) {
			if (p == null) continue;
			if (p.getSymbId() != null) ids.add(p.getSymbId());
		}
		
		return ids;
	}
	
	/**
	 * Get the source TTE data structure associated with a discrete model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.TimeToEventData
	 */
	public TimeToEventData getTimeToEventData() { return tte_data; }
	
	/**
	 * Flag if a discrete model has categorical data bound to it.
	 * @return boolean
	 * @see eu.ddmore.libpharmml.dom.modeldefn.CategoricalData
	 * 
	 */
	public boolean hasCategoricalData() { return categorical_data != null; }

	/**
	 * Report that a discrete model has a count data model.
	 * @return boolean
	 */
	public boolean hasCountData() { return count_data != null; }
	
	/**
	 * Flag if the observation model has local scoped parameters.
	 * @return boolean
	 */
	public boolean hasParameters() { return !param_refs.isEmpty(); }
	
	@Override
	public boolean hasSymbolId(String name) {
		boolean has = false;
		if (name != null) {
			for (CommonParameter s : params) {
				String v_name = s.getSymbId();
				if (v_name != null) {
					if (v_name.equals(name)) {
						has = true;
						break;
					}
				}
			}
		}
		
		return has;
	}
	
	/**
	 * Flag if the observation block has a TTE structure.
	 * @return boolean
	 * @see ObservationBlockImpl#getTimeToEventData()
	 */
	public boolean hasTimeToEventData() { return tte_data != null; }
	
	private void initErrorModel() {
		ContinuousObservationModel com = ob.getContinuousData();
		discrete = ob.getDiscrete();
		
		if (com != null) {
			List<PharmMLElement> elements = com.getListOfObservationModelElement();
			if (params != null) {
				for (PharmMLElement o : elements) {
					if (o == null) continue;

					if (isPopulationParameter(o)) {
						PopulationParameter p = (PopulationParameter) o;
						checkAssignment(p);
						String symbolId = p.getSymbId();
						if (!isSymbolIdOkay(symbolId)) throw new IllegalStateException("The node symbolId (" + symbolId + ") is not okay.");
						params.add((PopulationParameter) p);
					} else if (isRandomVariable(o)) {
						ParameterRandomVariable rv = (ParameterRandomVariable) o;
						String symbolId = rv.getSymbId();
						if (!isSymbolIdOkay(symbolId)) throw new IllegalStateException("The node symbolId (" + symbolId + ") is not okay.");
						params.add(rv);
						random_variables.add((ParameterRandomVariable) o);
					} else if (isIndividualParameter(o)) {
						IndividualParameter ip = (IndividualParameter) o;
						String symbolId = ip.getSymbId();
						if (!isSymbolIdOkay(symbolId)) throw new IllegalStateException("The node symbolId (" + symbolId + ") is not okay.");
						indiv_params.add(ip);
					} 
					else if (isLocalVariable(o)) continue; // Do nothing, unsupported model form.
					else if (isCorrelation(o)) register((Correlation) o); 
					else
						throw new UnsupportedOperationException("Observation scope parameter type unsupported ("+ o + ")");
				}
			}

			error_model = com.getObservationError();
		} else if (discrete != null) {
			count_data = discrete.getCountData();
			categorical_data = discrete.getCategoricalData();
			tte_data = discrete.getTimeToEventData();
		} else
			throw new NullPointerException("The observation model (blkId='" + ob.getBlkId() + "') has no bound observation error model.");
	}
	
	/**
	 * Test if the error model should be applied to a structural block.
	 * @param sb A Structural Block
	 * @return boolean
	 */
	public boolean isApplicable(StructuralBlock sb) {
		if (sb != null) {
			String sb_name = sb.getName();
			if (sb_name != null) {
				if (model_block_ids.contains(sb_name)) return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Flag if the error model is discrete.
	 * @return boolean
	 */
	public boolean isDiscrete() { return discrete != null; }
	
	/**
	 * Test if a parameter has error model (observation) scope.
	 * @param p Simple Parameter
	 * @return boolean
	 */
	public boolean isObservationParameter(PopulationParameter p) {
		if (p != null) {
			for (ObservationParameter param_ref : param_refs) {
				if (param_ref == null) continue;
				if (param_ref.param != null) {		
					if (param_ref.param.equals(p)) return true;
				}
			}
		}
		
		return false;
	}
	
	/**
	 * Check whether a model element is simulation output associated with an error model.
	 * @param v Model element under consideration.
	 * @return boolean
	 */
	public boolean isSimulationOutput(PharmMLRootType v) {
		if (v != null) {
			for (SimulationOutput simulation_output : simulation_outputs) {
				if (simulation_outputs == null) continue;
				if (simulation_output.v == null) continue;
				if (v.equals(simulation_output.v)) {
					return true;
				}
			}
		}
		
		return false;
	}

	/**
	 * Check whether the error model references a local variable in a calculation.
	 * @param v Local variable 
	 * @return boolean
	 */
	public boolean isUsing(VariableDefinition v) {
		if (v == null) return false;
		
		for (SimulationOutput output : simulation_outputs) {
			if (output == null) continue;
			if (v.equals(output.v)) return true;
		}
		
		return false;
	}
	
	@Override
	public String toString() { return getName(); }
}

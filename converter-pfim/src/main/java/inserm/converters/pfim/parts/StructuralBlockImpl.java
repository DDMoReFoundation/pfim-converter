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

import static crx.converter.engine.PharmMLTypeChecker.isDerivative;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.engine.VariableDeclarationContext;
import crx.converter.engine.assoc.Cluster;
import crx.converter.engine.common.DerivativeEvent;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.OrderableBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.CommonVariableDefinition;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.InitialCondition;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLElement;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.StandardAssignable;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modeldefn.StructuralModel;
import eu.ddmore.libpharmml.dom.modeldefn.pkmacro.PKMacro;
import eu.ddmore.libpharmml.pkmacro.translation.MacroOutput;

/**
 * Interpreted PharmMl structural block.
 */
public class StructuralBlockImpl extends PartImpl implements StructuralBlock, OrderableBlock {
	private static InitialCondition initialCondition(Double value) {
		InitialCondition ic = new InitialCondition();
		
		ic.setInitialValue(standardAssignable(value));
		ic.setInitialTime(standardAssignable(0.0));
		
		return ic;
	}
	
	private static StandardAssignable standardAssignable(double value) {
		StandardAssignable v = new StandardAssignable();
		v.assign(new RealValue(value));
		return v;
	}
	
	private Accessor a = null;
	private List<BinaryTree> bts =  new ArrayList<BinaryTree>();
	private List<PharmMLRootType> cached_declaration_list = new ArrayList<PharmMLRootType>(); 
	private VariableDefinition doseTimingVariable = null;
	private List<DerivativeEvent> events = new ArrayList<DerivativeEvent>();
	private Map<String, VariableDefinition> local_map_name = new HashMap<String, VariableDefinition>();
	private List<VariableDefinition> locals = new ArrayList<VariableDefinition>();
	private List<PKMacro> macro_list = null;
	private MacroOutput macro_output = null;
	private List<PopulationParameter> params = new ArrayList<PopulationParameter>();
	private StructuralModel sm = null;
	private Map<CommonVariableDefinition, Integer> state_map_idx = new HashMap<CommonVariableDefinition, Integer>();
	private Map<String, CommonVariableDefinition> state_map_name = new HashMap<String, CommonVariableDefinition>();
	private List<DerivativeVariable> states = new ArrayList<DerivativeVariable>();
	private TreeMaker tm = null;
	private Map<Object, BinaryTree> tree_map = new HashMap<Object, BinaryTree>();
	
	/**
	 * Constructor
	 * @param sm_ Structural Model
	 * @param lexer_ Converter Instance
	 */
	public StructuralBlockImpl(StructuralModel sm_, ILexer lexer_) {
		if (sm_ == null) throw new NullPointerException("The structural model cannot be NULL.");
		sm = sm_;
		if (sm.getBlkId() == null) throw new IllegalStateException("State model must have a BLK ID");
		
		lexer = lexer_;
		tm = lexer.getTreeMaker();
		a = lexer.getAccessor();
		
		Integer idx = 0;
		for (Object v : sm.getListOfStructuralModelElements()) {
			if (isDerivative(v)) {
				DerivativeVariable dv = (DerivativeVariable) v;
				String symbolId = dv.getSymbId();
				if (!isSymbolIdOkay(symbolId)) throw new IllegalStateException("The node symbolId (" + symbolId + ") is not okay.");
				idx++;
				
				if (lexer.isUsePiecewiseAsEvents()) {
					// Infers a conditional statement statement on a state variable
					if (hasPiecewise(dv)) {
						DerivativeEvent evt = new DerivativeEvent(dv, a);
						events.add(evt);
					}
				}
				
				state_map_idx.put(dv, idx);
				state_map_name.put(symbolId, dv);
				states.add((DerivativeVariable) dv);
			}
		}
		
		for (Object v : sm.getListOfStructuralModelElements()) {
			if (isLocalVariable(v)) {
				VariableDefinition local = (VariableDefinition) v;
				String symbolId = local.getSymbId();
				if (!isSymbolIdOkay(symbolId)) throw new IllegalStateException("The node symbolId (" + symbolId + ") is not okay.");
				
				// Check if a local variable is a dose variable specified by an external file reference.
				// Isolate as required from the locals list.
				boolean isConditionalDoseEventTarget = false;
				
				VariableDeclarationContext ctx = lexer.guessContext(local);
				
				if (lexer.isIsolatingDoseTimingVariable()) {
					if (VariableDeclarationContext.DT.equals(ctx)) {
						doseTimingVariable = local;
						continue;
					} else if (VariableDeclarationContext.DOSE.equals(ctx)) continue;
				} 
				
				if (isConditionalDoseEventTarget && lexer.isIsolateConditionalDoseVariable()) continue; // Skip so it doesn't get added to the locals list.
				
				locals.add(local);
				local_map_name.put(symbolId, local);
			}
		}
		
		for (Object v : sm.getListOfStructuralModelElements()) if (isPopulationParameter(v)) params.add((PopulationParameter) v);
	}
	
	private void add(Object key, BinaryTree value) {
		if (key != null && value != null) {
			if (!tree_map.containsKey(key)) tree_map.put(key, value);
		}
	}
	
	@Override
	public boolean addCluster(Cluster cluster) { throw new UnsupportedOperationException(); }
	
	@Override
	public void buildTrees() {
		if (treeBuilt) return;
		
		bts.clear();
		
		createStateTrees();
		createLocalVariableTrees();
		createInitialConditionTrees();
		createParameterTrees();
		treeBuilt = true;
	}
	
	
	private void clearAllVariableLists() {
		params.clear();
		locals.clear();
		states.clear();
	}
	
	/**
	 * Check that the structural model has the model element.
	 * @param v
	 * @return boolean
	 */
	public boolean contains(PharmMLRootType v) {
		if (isDerivative(v)) {
			if (states.contains(v)) return true;
		} else if (isLocalVariable(v)) {
			if (locals.contains(v)) return true;
		}
		
		return false;
	}
	
	private void createInitialConditionTrees() {
		for (DerivativeVariable s : states) {
			if (s.getHistory() != null) continue; // DDE model so no initial condition
			
			if (s.getInitialCondition() == null) s.setInitialCondition(initialCondition(0.0));
			InitialCondition ic = s.getInitialCondition();
			
			BinaryTree bt = tm.newInstance(ic);
			bts.add(bt);
			lexer.addStatement(ic, bt);
			lexer.updateNestedTrees();
		}
	}
	
	private void createLocalVariableTrees() {
		for (VariableDefinition v : locals) {
			BinaryTree bt = tm.newInstance(v);
			lexer.updateNestedTrees();
			bts.add(bt);
			add(v, bt);
			lexer.addStatement(v, bt);
		}
	}
	
	private void createParameterTrees() {
		if (params.isEmpty()) return;
		
		for (PopulationParameter p : params) {
			if (p == null) continue;
			
			BinaryTree bt = tm.newInstance(p);
			lexer.addStatement(p, bt);
		}
	}
	
	private void createStateTrees() {
		for (DerivativeVariable s : states) {
			BinaryTree bt = tm.newInstance(s);
			lexer.updateNestedTrees();
			bts.add(bt);
			add(s, bt);
			lexer.addStatement(s, bt);
		}
	}
	
	@Override
	public List<Cluster> getClusters() { throw new UnsupportedOperationException(); }
	
	/**
	 * Get the dose timing (DT) variable if defined in the model.
	 * @return VariableDefinition
	 */
	public VariableDefinition getDoseTimingVariable() { return doseTimingVariable; }

	/**
	 * Get a list of structural model event.
	 * @return java.util.List<DerivativeEvent>
	 */
	public List<DerivativeEvent> getEvents() { return events; }
	
	/**
	 * All of the declared variables in the structural model.
	 * @return List<PharmMLRootType>
	 */
	public List<PharmMLRootType> getListOfDeclarations() {
		if (!cached_declaration_list.isEmpty()) return cached_declaration_list;
		
		List<PharmMLRootType> variables = new ArrayList<PharmMLRootType>();
	
		variables.addAll(getParameters());
		variables.addAll(getLocalVariables());
		variables.addAll(getStateVariables());
		
		return variables;
	}
	
	/**
	 * List of local variables
	 * @return java.util.List<eu.ddmore.libpharmml.dom.commontypes.VariableDefinition>
	 */
	public List<VariableDefinition> getLocalVariables() { return locals; }
	
	/**
	 * Get the source structural model.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.StructuralModel
	 */
	public StructuralModel getModel() { return sm; }
	
	@Override
	public String getName() {
		String blkId = "_RUBBISH_DEFAULT_VALUE_";
		if (sm != null) {
			if (sm.getBlkId() != null) blkId = sm.getBlkId();
		}
		return blkId;
	}
		
	/**
	 * Get the list of structural parameters.
	 * @return java.util.List<eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter>
	 */
	public List<PopulationParameter> getParameters() { return params; }
	
	/**
	 * Get the raw macro outputs generated by the translator component.
	 * @see eu.ddmore.libpharmml.pkmacro.translation.Translator#translate
	 * @return eu.ddmore.libpharmml.pkmacro.translation.MacroOutput
	 */
	public MacroOutput getPKMacroOutput() { return macro_output ; }
	
	/**
	 * Get the source macros of the structural model.
	 * @return List<eu.ddmore.libpharmml.dom.modeldefn.pkmacro.PKMacro>
	 */
	public List<PKMacro> getPKMacros() { return macro_list; }
	
	@Override
	public Integer getStateVariableIndex(String name) {
		Integer idx = -1;
		if (name == null) return idx;
		
		if (state_map_name.containsKey(name)) {
			CommonVariableDefinition s = (state_map_name.get(name));
			if (s != null) {
				if (state_map_idx.containsKey(s)) idx = state_map_idx.get(s);
			}
		}

		return idx; 
	}
	
	/**
	 * Get a list of derivatvies
	 * @return java.util.List<eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable>
	 */
	public List<DerivativeVariable> getStateVariables() { return states; }
	
	@Override
	public List<String> getSymbolIds() {
		ArrayList<String> ids = new ArrayList<String>();
		for (DerivativeVariable s : states) {
			if (s == null) continue;
			if (s.getSymbId() != null) ids.add(s.getSymbId());
		}
		
		return ids;
	}
	
	/**
	 * Flag if a dose timing (DT) variable is defined in the model
	 * @return boolean
	 */
	public boolean hasDoseTimingVariable() { return doseTimingVariable != null; }
	
	/**
	 * Flag that the structural model has simulation events.
	 * @return boolean
	 */
	public boolean hasEvents() { return events.size() > 0; }
	
	/**
	 * Flag that the structural model has parameters.
	 * @return boolean
	 */
	public boolean hasParameters() { return !params.isEmpty(); }
	
	private boolean hasPiecewise(DerivativeVariable dv) {
		if (dv != null) {
			Rhs eq = dv.getAssign();
			if (eq != null) if (eq.getPiecewise() != null) return true;
		}
		
		return false;
	}
	
	/**
	 * Flag if the structural model has source PK macros.
	 * @return boolean
	 */
	public boolean hasPKMacros() {
		if (macro_list == null) return false;
		else return macro_list.size() > 0;
	}
	
	@Override
	public boolean hasSymbolId(String name) {
		if (name != null) {
			for (DerivativeVariable s : states) {
				String v_name = s.getSymbId();
				if (v_name == null) continue;
				else if (v_name.equals(name)) return true;
			}
			
			for (VariableDefinition local : locals) {
				String v_name = local.getSymbId();
				if (v_name == null) continue;
				else if (v_name.equals(name)) return true;
			}
			
			for (PopulationParameter p : params) {
				String v_name = p.getSymbId();
				if (v_name == null) continue;
				else if (v_name.equals(name)) return true;
			}	 
		}
		
		return false;
	}
	
	/**
	 * Flag that the model is a delayed-event ODE.
	 * @return boolean
	 */
	public boolean isDDE() { return false; }
	
	@Override
	public boolean isMixedEffect() { throw new UnsupportedOperationException(); }
	
	/**
	 * Flag that model is a standard ODE form.
	 * @return boolean
	 */
	public boolean isODE() { return !states.isEmpty(); }
	
	/**
	 * Flag that model is plain function, i.e. no derivatives.
	 * @return boolean
	 */
	public boolean isPlainFunction() { return states.size() == 0; }
	
	/**
	 * Test if variable is a derivative.
	 * @param name Variable name
	 * @return boolean
	 */
	public boolean isStateVariable(String name) {
		boolean dt_something = false;
		
		if (name != null) {
			for (DerivativeVariable state : states) {
				if (state != null) {
					String symId = state.getSymbId();
					if (symId != null) {
						if (symId.equals(name)) {
							dt_something = true;
							break;
						}
					}
				}
			}
		}
		
		return dt_something;
	}
	
	/**
	 * Flag if the structural block is using untranslated PK Macros.
	 * @return boolean
	 */
	public boolean isUsingUntranslatedPKMacros() { throw new UnsupportedOperationException(); }
	
	/**
	 * Set the ordered variable list within the structural block.<br/>
	 * This is set outside of the StructuralBlock, hence this accessor function.
	 * @param ordered_variables Ordered variable List.
	 */
	public void setOrderedVariableList(List<PharmMLElement> ordered_variables) {
		if (ordered_variables == null) return;
		if (ordered_variables.isEmpty()) return;
		
		List<PopulationParameter> old_params = new ArrayList<PopulationParameter>();
		List<VariableDefinition> old_locals = new ArrayList<VariableDefinition>();
		List<DerivativeVariable> old_states = new ArrayList<DerivativeVariable>(); 
		
		old_params.addAll(params);
		old_locals.addAll(locals);
		old_states.addAll(states);
		
		clearAllVariableLists();
		
		for (PharmMLElement ordered_variable : ordered_variables) {
			if (isPopulationParameter(ordered_variable)) {
				PopulationParameter p = (PopulationParameter) ordered_variable;
				if (old_params.contains(p)) params.add(p);
			} else if (isLocalVariable(ordered_variable)) {
				VariableDefinition v = (VariableDefinition) ordered_variable;
				if (old_locals.contains(v)) locals.add(v);
			} else if (isDerivative(ordered_variable)) {
				DerivativeVariable dv = (DerivativeVariable) ordered_variable;;
				if (old_states.contains(dv)) states.add(dv);
			}
		}
		
		
		cached_declaration_list.clear();
		for (PharmMLElement o : ordered_variables) {
			if (isPopulationParameter(o) || isLocalVariable(o) || isDerivative(o)) cached_declaration_list.add((PharmMLRootType) o);
		}
	} 
	
	/**
	 * Get the raw macro outputs generated by the translator component.
	 * @param macro_output_
	 * @see eu.ddmore.libpharmml.pkmacro.translation.Translator#translate
	 */
	public void setPKMacroOutput(MacroOutput macro_output_) { macro_output = macro_output_; }

	/**
	 * Set the PK macro list, as read from a parent structural model.
	 * @param macro_list_ PK Macro list
	 */
	public void setPKMacros(List<PKMacro> macro_list_) { macro_list = macro_list_; }

	@Override
	public List<CommonVariableDefinition> getRegressors() { throw new UnsupportedOperationException(); }

	@Override
	public boolean isRegressor(CommonVariableDefinition v) { throw new UnsupportedOperationException(); }
}

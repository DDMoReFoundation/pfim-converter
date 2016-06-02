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

import static crx.converter.engine.PharmMLTypeChecker.isVariabilityLevelDefinition;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.spi.ILexer;
import crx.converter.spi.blocks.VariabilityBlock;
import eu.ddmore.libpharmml.dom.commontypes.LevelReference;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.ParentLevel;
import eu.ddmore.libpharmml.dom.modeldefn.Variability;
import eu.ddmore.libpharmml.dom.modeldefn.VariabilityDefnBlock;
import eu.ddmore.libpharmml.dom.modeldefn.VariabilityLevelDefinition;

/**
 * Wrapper class for a PharmML variability block.
 */
public class VariabilityBlockImpl extends PartImpl implements VariabilityBlock  {
	private Map<String, StringValue> comparator_map = new HashMap<String, StringValue>();
	private boolean parameter_variability = false;
	private boolean residual_error = false;
	private Map<String, String> level_dependency_map = new HashMap<String, String>();
	private List<String> symbols = new ArrayList<String>();
	private VariabilityDefnBlock vb = null;
	
	/**
	 * Constructor 
	 * @param vb_ Source Model
	 */
	public VariabilityBlockImpl(VariabilityDefnBlock vb_, ILexer lexer_) {
		if (vb_ == null) throw new NullPointerException("VariabilityDefnBlock argument cannot be null.");
		if (lexer_ == null) throw new NullPointerException("Converter/Lexer argument cannot be null.");
		
		vb = vb_;
		lexer = lexer_;
		Variability type = vb.getType();
		
		if (type == Variability.PARAMETER_VARIABILITY) parameter_variability = true;
		else if (type == Variability.RESIDUAL_ERROR) residual_error = true;
		else throw new UnsupportedOperationException("Variability not supported (" + type + ")");
		
		Accessor a = lexer.getAccessor();
		List<VariabilityLevelDefinition> levels = vb.getLevel();
		
		for (VariabilityLevelDefinition level : levels) {
			String symbId = level.getSymbId();
			if (symbId == null) throw new NullPointerException("Level ID is NULL.");
			else symbols.add(symbId);
			
			if (level.getParentLevel() != null) {
				ParentLevel pl = level.getParentLevel();
				SymbolRef ref = pl.getSymbRef();
				PharmMLRootType element = a.fetchElement(pl.getSymbRef());
				if (!isVariabilityLevelDefinition(element)) 
					throw new IllegalStateException("The parent model element is not a variability level (symbId='" + ref.getSymbIdRef() + "')");
				
				if (level_dependency_map.containsKey(symbId)) throw new IllegalStateException("A variability level appears to have cyclic referencing.");
				level_dependency_map.put(symbId, ref.getSymbIdRef());
			}
			
			addToCompareMap(symbId);
		}
	}  
	
	private void addToCompareMap(String symbol) {
		if (symbol == null) return;
		if (!comparator_map.containsKey(symbol))  comparator_map.put(symbol, new StringValue(symbol));
	}
	
	@Override
	public void buildTrees() { }
	
	/**
	 * See if a putative child level is support by a parent level.
	 * @param levelChild Child Level
	 * @param levelParent Parent Level
	 * @return boolean
	 */
	public boolean dependsUpon(String levelChild, String levelParent) {
		if (levelChild == null || levelParent == null) return false;
		if (!hasSymbolId(levelChild) || !hasSymbolId(levelParent)) return false;
		
		return false;
	}
	
	/**
	 * Get the levels declared in the variability scope.
	 * @return List<String>
	 */
	public List<String> getLevels() { return symbols; }
	
	/**
	 * Get the source model element.
	 * @return eu.ddmore.libpharmml.dom.modeldefn.VariabilityDefnBlock
	 */
	public VariabilityDefnBlock getModel() { return vb; }
	
	@Override
	public String getName() {
		if (vb != null) return vb.getBlkId();
		else return null;
	}
	
	@Override
	/**
	 * In this case the symbol IDs of the block are the variability levels.
	 * @return List<String>
	 */
	public List<String> getSymbolIds() { return symbols; }

	/**
	 * Check if a random variable has scope in the variability block.
	 * @param rv Random variable
	 * @return boolean
	 */
	public boolean hasScope(ParameterRandomVariable rv) {
		if (rv == null) return false;
		
		List<LevelReference> lrefs = rv.getListOfVariabilityReference();
		
		if (lrefs == null) return false;
		for (LevelReference lref : lrefs) {
			if (lref == null) continue;
			SymbolRef ref = lref.getSymbRef();
			if (ref != null) {
				if (hasSymbolId(ref.getSymbIdRef())) return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Check if a random variable has scope in the variability block for a specific level.
	 * @param level Level Name
	 * @param rv Random variable
	 * @return boolean
	 */
	public boolean hasScope(String level, ParameterRandomVariable rv) {
		if (rv == null || level == null) return false;
		if (!hasSymbolId(level)) return false; 
		
		List<LevelReference> lrefs = rv.getListOfVariabilityReference();
		for (LevelReference lref : lrefs) {
			if (lref == null) continue;
			
			SymbolRef ref = lref.getSymbRef();
			if (ref == null) continue;
			
			String referenced_level = ref.getSymbIdRef();
			if (hasSymbolId(referenced_level)) return true;
		}
		
		return false;
	}
	
	@Override
	public boolean hasSymbolId(String name) {
		if (name != null) return symbols.contains(name);
		return false;
	}
	
	/**
	 * Flag if block linked to parameter variability scope, i.e. individual random things.
	 * @return boolean
	 */
	public boolean isParameterVariability() { return parameter_variability; }
	

	/**
	 * Flag if block linked to residual scope, i.e. linked to time series error model.
	 * @return boolean
	 */
	public boolean isResidualError() { return residual_error; }
	
	/**
	 * Read the level/scoping string bound to a random variable.
	 * @param rv Random Variable
	 * @return String Level String
	 */
	public String readLevel(ParameterRandomVariable rv) {
		String level = null;
		
		if (rv != null) {
			List<LevelReference> lrefs = rv.getListOfVariabilityReference();
			if (lrefs != null) {
				if (lrefs.size() > 0) {
					LevelReference lref = lrefs.get(0);
					if (lref != null) {
						SymbolRef ref = lref.getSymbRef();
						if (ref != null) level = ref.getSymbIdRef();
					}
				}
			}
		}
		
		return level;
	}
	
	@Override
	public String toString() { return getName(); }
}

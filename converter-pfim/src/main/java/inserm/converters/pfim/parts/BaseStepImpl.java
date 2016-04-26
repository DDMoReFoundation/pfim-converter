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

import static crx.converter.engine.PharmMLTypeChecker.isColumnMapping;
import static crx.converter.engine.PharmMLTypeChecker.isCovariate;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isMultipleDVMapping;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.engine.CategoryRef_;
import crx.converter.engine.common.ConditionalDoseEvent;
import crx.converter.engine.common.ConditionalDoseEventRef;
import crx.converter.engine.common.MultipleDvRef;
import crx.converter.engine.common.TemporalDoseEvent;
import crx.converter.engine.common.BaseTabularDataset.ElementMapping;
import crx.converter.spi.steps.BaseStep;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.SymbolType;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.dataset.CategoryMapping;
import eu.ddmore.libpharmml.dom.dataset.ColumnDefinition;
import eu.ddmore.libpharmml.dom.dataset.ColumnMapping;
import eu.ddmore.libpharmml.dom.dataset.ColumnReference;
import eu.ddmore.libpharmml.dom.dataset.ColumnType;
import eu.ddmore.libpharmml.dom.dataset.DataSet;
import eu.ddmore.libpharmml.dom.dataset.HeaderColumnsDefinition;
import eu.ddmore.libpharmml.dom.dataset.IgnoreLine;
import eu.ddmore.libpharmml.dom.dataset.MapType;
import eu.ddmore.libpharmml.dom.maths.Piece;
import eu.ddmore.libpharmml.dom.maths.Piecewise;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;
import eu.ddmore.libpharmml.dom.trialdesign.ExternalDataSet;
import eu.ddmore.libpharmml.dom.trialdesign.MultipleDVMapping;


/**
 * Abstract/Base analysis step.
 * Methods mostly to process an external data set references and column mappings.
 */
public abstract class BaseStepImpl extends PartImpl implements BaseStep {
	/**
	 * Accessor instance.
	 */
	protected Accessor a = null;
	
	//private String categorical_dose_event_colName = null;
	private List<ConditionalDoseEvent> category_dose_evts = new ArrayList<ConditionalDoseEvent>();
	
	/**
	 * Column name map.
	 */
	protected Map<String, ColumnDefinition> column_name_map = null;
	/**
	 * List of conditional dose events.
	 */
	protected List<ConditionalDoseEvent> conditional_dose_events = new ArrayList<ConditionalDoseEvent>();
	
	private List<VariableDefinition> conditional_dose_evt_targets = new ArrayList<VariableDefinition>();
	private Map<VariableDefinition, ConditionalDoseEventRef> dose_ref_var_map = new HashMap<VariableDefinition, ConditionalDoseEventRef>();
	
	/**
	 * Element mappings list.
	 */
	protected List<ElementMapping> element_mappings = new ArrayList<ElementMapping>();
	
	/**
	 * External dataset.
	 */
	protected ExternalDataSet exd = null;
	
	private String ignoreLineSymbol = null;
	private ColumnDefinition infusion_column = null;
	private List<MultipleDvRef> mdv_refs = new ArrayList<MultipleDvRef>();
	private TemporalDoseEvent tde = null;
	
	/**
	 * Add an element mapping from an external data file.
	 * @param mapping Element Mapping
	 */
	protected void addElementMapping(ElementMapping mapping) { if (mapping != null) element_mappings.add(mapping); }
	
	/**
	 * Build the data set mappings bound to the current analysis step.
	 */
	protected void buildExternalDatasetMappings() {
		if (exd == null) return;
		
		// Follow the newer external dataset way of doing things.
		DataSet ds = exd.getDataSet();
		createColumnNameMap(ds);
		
		// Populate the element mappings of the dataset to the model.
		for (PharmMLRootType o :  exd.getListOfColumnMappingOrColumnTransformationOrMultipleDVMapping()) {
			if (isColumnMapping(o)) processColumnMapping((ColumnMapping) o);
			else if (isMultipleDVMapping(o)) processMDVMapping((MultipleDVMapping) o);
			else throw new UnsupportedOperationException("Data column bound with an unknown context (not support yet, object='" + o +  "')");
		}
	}
	
	/**
	 * Create an internal map of column definition to column name.
	 * @param ds Dataset containing the column definitions.
	 */
	protected void createColumnNameMap(DataSet ds) {
		if (ds == null) throw new NullPointerException("The dataset is NULL");
		
		column_name_map = new HashMap<String, ColumnDefinition>();
		
		HeaderColumnsDefinition def = ds.getDefinition();
		if (def == null) throw new NullPointerException("Column Definition list is NULL");
		List<ColumnDefinition> cols = def.getListOfColumn();
		if (cols == null) throw new NullPointerException("No columns defined.");
		if (cols.isEmpty()) throw new IllegalStateException("No data columns defined in the objective data set.");
		
		for (ColumnDefinition col : cols) {
			if (col == null) continue;
			String name = col.getColumnId();
			Integer col_idx = col.getColumnNum();
			ColumnType usage = col.getListOfColumnType().get(0);
			SymbolType datatype = col.getValueType();
			
			if (name == null) throw new NullPointerException("Dataset column name is NULL");
			if (col_idx == null) throw new NullPointerException("Dataset column index is unspecified");
			
			// Only check data type is column is being used for something, otherwise ignore as not used by software.
			if (usage != null) {
				if (datatype == null) throw new IllegalStateException("The column data type is not specified");
				if (ColumnType.RATE.equals(usage)) {
					infusion_column = col;
				}
			}
			
			if (column_name_map.containsKey(name)) throw new IllegalStateException("Column name map contains duplicate key.");
			column_name_map.put(name, col);
		}
	}
	
	/**
	 * Get the list of categorical dose events as mapped to a single column in 
	 * an external dataset.<br/>
	 * Basically this flags the NONMEM/Monolix dosing target in a data frame.
	 * @return List<ConditionalDoseEvent>
	 */
	public List<ConditionalDoseEvent> getCategoricalDoseEvents() { return category_dose_evts; }
	
	/**
	 * Get a named column associated with an external data set.
	 * @param name Column Name
	 * @return eu.ddmore.libpharmml.dom.dataset.ColumnDefinition
	 */
	public ColumnDefinition getColumn(String name) {
		if (name == null) return null;
		if (column_name_map.containsKey(name)) return column_name_map.get(name); 
		else return null;
	}
	
	/**
	 * Get the columns for a specific usage.
	 * Returns a list as dose possible for multiple compartments or multiple covariate columns.
	 * @param usage Usage
	 * @return java.util.List<ElementMapping>
	 */
	public List<ColumnDefinition> getColumns(ColumnType usage) {
		List<ColumnDefinition> list = new ArrayList<ColumnDefinition>();
		
		if (usage != null) {
			for (ColumnDefinition col : column_name_map.values()) {
				if (col == null) continue;
				ColumnType type = col.getListOfColumnType().get(0);
				if (type.equals(usage)) {
					if (!list.contains(col)) list.add(col);
				}
			}
		}
		
		return list;
	}
	
	/**
	 * Get the conditional dose event associated with a variable declaration.
	 * @return ConditionalDoseEventRef
	 */
	public ConditionalDoseEventRef getConditionalDoseEventRef(VariableDefinition v) {
		if (v != null) if (dose_ref_var_map.containsKey(v)) return dose_ref_var_map.get(v);
		return null;
	}
	
	/**
	 * Get a list of conditional dose events as read from an external data file.
	 * @return java.util.List<ConditionalDoseEventRef>
	 */
	public List<ConditionalDoseEventRef> getConditionalDoseEventRefs() {
		List<ConditionalDoseEventRef> list = new ArrayList<ConditionalDoseEventRef>();
		
		for (ConditionalDoseEventRef ref : dose_ref_var_map.values()) if (!list.contains(ref)) list.add(ref);
		
		return list;
	}
	
	/**
	 * Get the conditional dose events (if any) linked to the estimation.
	 * @return java.util.List<ConditionalDoseEvent>
	 */
	public List<ConditionalDoseEvent> getConditionalDoseEvents() { return conditional_dose_events; }
	
	/**
	 * Get the element mapping to the objective dataset.
	 * @return java.util.List<ElementMapping>
	 */
	public List<ElementMapping> getElementMappings() { return element_mappings; }
	
	/**
	 * Get the element mapping to the objective dataset.
	 * @param usage Column Usage
	 * @return java.util.List<ElementMapping>
	 */
	public List<ElementMapping> getElementMappings(ColumnType usage) {
		List<ElementMapping> list = new ArrayList<ElementMapping>();
		
		if (usage != null) {
			for (ElementMapping element_mapping : element_mappings) {
				if (element_mapping == null) continue;
				ColumnDefinition col = element_mapping.getColumnDefinition();
				if (col == null) continue;
				if (usage.equals(col.getListOfColumnType().get(0))) list.add(element_mapping);
			}
		}
		
		return list; 
	}
	
	/**
	 * Return specific model elements associated with a column usage.
	 * @param usage
	 * @return java.util.List<eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType>
	 */
	public List<PharmMLRootType> getElements(ColumnType usage) {
		List<PharmMLRootType> list = new ArrayList<PharmMLRootType>();
		
		if (usage != null) {
			for (ElementMapping element_mapping : element_mappings) {
				if (element_mapping == null) continue;
				ColumnDefinition col = element_mapping.getColumnDefinition();
				if (col != null) {
					ColumnType type = col.getListOfColumnType().get(0);
					if (usage.equals(type)) {
						PharmMLRootType element = element_mapping.getElement();
						if (element != null && !list.contains(element)) list.add(element);
					}
				}
			}
		}
		
		return list;
	}
	
	/**
	 * Get the External data set that contains the objectivve data file reference.
	 * @return eu.ddmore.libpharmml.dom.modellingsteps.ExternalDataSet
	 */
	public ExternalDataSet getExternalDataSet () {
		return exd;
	}
	
	/**
	 * Get the ignore line/row character symbol for a row in a data frame.<br/>
	 * @return java.lang.String.
	 */
	public String getIgnoreLineSymbol() {
		if (ignoreLineSymbol == null) {
			if (exd != null) {
				DataSet ds = exd.getDataSet();
				if (ds != null) {
					HeaderColumnsDefinition definition  = ds.getDefinition();
					if (definition != null) {
						IgnoreLine igl = definition.getIgnoreLine();
						if (igl != null) {
							ignoreLineSymbol = igl.getSymbol();
						}
					}
				}
			}
		}
		
		return ignoreLineSymbol;
	}
	
	/**
	 * Get the column definition for the declared infusion column.
	 * @return ColumnDefinition
	 */
	public ColumnDefinition getInfusionColumn() { return infusion_column;}
	
	/**
	 * Flag that the step has a multiple DV mappings references mapped to an external data file.
	 * This construct can be used to switch a DV column to different error models based on the value read from a CSV file.
	 * @return boolean
	 * @see eu.ddmore.libpharmml.dom.trialdesign.MultipleDVMapping
	 * @see eu.ddmore.libpharmml.dom.trialdesign.ExternalDataSet#getListOfColumnMappingOrColumnTransformationOrMultipleDVMapping() 
	 */
	public List<MultipleDvRef> getMultipleDvRefs() { return mdv_refs; }
	
	/**
	 * Get Temporal Dose event as read from an external dataset.
	 * @return TemporalDoseEvent
	 */
	public TemporalDoseEvent getTemporalDoseEvent() { return tde; }
	
	/**
	 * Flag if estimation involves conditional dosing.
	 * @return boolean
	 * @see ConditionalDoseEvent
	 */
	public boolean hasConditionalDoseEvents() { return conditional_dose_events.size() > 0; }
	
	/**
	 * Flag that the step has a multiple DV mappings references mapped to an external data file.
	 * This construct can be used to switch a DV column to different error models based on the value read from a CSV file.
	 * @return boolean
	 * @see eu.ddmore.libpharmml.dom.trialdesign.MultipleDVMapping
	 * @see eu.ddmore.libpharmml.dom.trialdesign.ExternalDataSet#getListOfColumnMappingOrColumnTransformationOrMultipleDVMapping() 
	 */
	public boolean hasMultipleDVRefs() { return mdv_refs.size() > 0; }
	
	/**
	 * Flag whether the modelling the Temporal Dose event as read from an external dataset.
	 * @return boolean
	 */
	public boolean hasTemporalDoseEvent() { return tde != null; }
	
	@Override
	public boolean isCategoricalDoseTargetColumn(String colName) { throw new UnsupportedOperationException(); }
	
	/**
	 * Check if a local variable is a conditional dosing target variable.<br/>
	 * If so, then needs to be isolated in languages like NONMEM.
	 * @param v
	 * @return boolean
	 */
	public boolean isConditionalDoseEventTarget(VariableDefinition v) {
		if (v == null) return false;
		else return conditional_dose_evt_targets.contains(v);
	}
	
	/**
	 * Flag is a step/task is undergoing an infusion process
	 * @return boolean
	 */
	public boolean isDoingInfusion() { throw new UnsupportedOperationException(); }
	
	/**
	 * Process a column mapping of a model element to a input column.
	 * @param mapping Column Mapping
	 */
	protected void processColumnMapping(ColumnMapping mapping) {
		ColumnReference colref = mapping.getColumnRef();
		if (colref == null) throw new NullPointerException("Mapping column reference is NULL");
		
		String colName = colref.getColumnIdRef();
		if (!column_name_map.containsKey(colName)) throw new IllegalStateException("Data set does not contain mapped data column (col='" + colName +"')");
		
		if (mapping.getSymbRef() != null) {
			PharmMLRootType element = a.fetchElement(mapping.getSymbRef());
			ElementMapping em = new ElementMapping(column_name_map.get(colName), element);
			addElementMapping(em);
			
			// Register category mappings (if any). 
			if (isCovariate(element)) registerCategoryMappings((CovariateDefinition) element, mapping);
			
			// Associate the category mappings to the elememt mapping record.
			if (!mapping.getListOfCategoryMapping().isEmpty()) {
				CategoryMapping cat_mappings = mapping.getListOfCategoryMapping().get(0);
				if (cat_mappings != null) {
					for (MapType cat_map : cat_mappings.getListOfMap()) em.addCategoryMappings(cat_map);
				}
			}
		}
		else if (mapping.getPiecewise() != null) {
			ColumnDefinition col = column_name_map.get(colName);
			// Check the context of usage.
			if (col != null) {
				ColumnType cxt = col.getListOfColumnType().get(0);
				if (cxt == null) throw new NullPointerException("The column context of usage is NULL.");
				
				if (cxt.equals(ColumnType.DOSE)) processDosePiecewiseMappings(col, mapping.getPiecewise());
				else if (cxt.equals(ColumnType.IDV) || cxt.equals(ColumnType.TIME)) processTemporalDosePiecewiseMappings(col, mapping.getPiecewise());
			}
		} 
	}
	
	private boolean processDosePiecewiseMappings(ColumnDefinition col, Piecewise pw) {
		List<Piece> pieces = pw.getListOfPiece();
		if (pieces.isEmpty()) throw new IllegalStateException("Piecewise block contains no conditional statements");
	
		Integer idx = 1;
		Accessor a = lexer.getAccessor();
		
		List<ConditionalDoseEvent> evts_copy = new ArrayList<ConditionalDoseEvent>();
		for (Piece piece : pieces) {
			ConditionalDoseEvent evt  = new ConditionalDoseEvent(lexer, col, piece);
			evts_copy.add(evt);
			
			PharmMLRootType element = evt.getTargetElement();
			if (element != null) {
				// Dosing on a local so add a reference so that a value can be 
				// assigned be externally assigned to the model symbol.
				if (isLocalVariable(element)) {
					VariableDefinition v = (VariableDefinition) element;
					conditional_dose_evt_targets.add(v);
					
					// Create a global reference to a dosing variable.
					if (lexer.isUseGlobalConditionalDoseVariable()) {
						String format = "evt_%s_%s_%s";
						String evt_id = String.format(format, col.getColumnId(), v.getSymbId(), idx);
					
						ConditionalDoseEventRef ref = new ConditionalDoseEventRef(evt, evt_id);
						a.register(evt_id, ref);
						v.assign(ref);
						dose_ref_var_map.put(v, ref);
					
						idx++;
					}
				}
			}
			
			conditional_dose_events.add(evt);
		}
		
		// If more than copy event, multiple pieces of is a category mapping for a dose event.
		if (evts_copy.size() > 1) category_dose_evts.addAll(evts_copy);
			
		return true;
	}
	
	private boolean processTemporalDosePiecewiseMappings(ColumnDefinition col, Piecewise pw) {
		List<Piece> pieces = pw.getListOfPiece();
		if (pieces.isEmpty()) throw new IllegalStateException("Piecewise block contains no conditional statements");
		
		Piece piece = pieces.get(0);
		if (piece == null) return false;
		
		ConditionalDoseEvent evt  = new ConditionalDoseEvent(lexer, col, piece);
		PharmMLRootType element = evt.getTargetElement();
		if (element == null) return false;
			
		if (isLocalVariable(element)) tde = new TemporalDoseEvent(lexer, col, piece);
		
		return true;
	}
	
	private void registerCategoryMappings(CovariateDefinition v, ColumnMapping mapping) {
		if (v == null || mapping == null) return;
		
		String blkId = a.getBlockId(v);
		
		List<CategoryMapping> cat_mappings = mapping.getListOfCategoryMapping();
		if (cat_mappings == null) return;
		if (cat_mappings.isEmpty()) return;
		
		List<CategoryRef_> element_category_list = new ArrayList<CategoryRef_>();
		for (CategoryMapping item : cat_mappings) {
			if (item == null) continue;
			for (MapType map : item.getListOfMap()) {
				if (map == null) continue;
				
				CategoryRef_ ref = new CategoryRef_(blkId, map.getModelSymbol(), map.getDataSymbol());
				ref.setElement(v);
				
				a.register(ref);
				element_category_list.add(ref);
			}
		}
		
		if (element_category_list.size() > 0) a.register(v, element_category_list);
	}
	
	/**
	 * Process the MDV mapping of an data input column to model elements,
	 * @param mdv MDV mapping definition
	 */
	protected void processMDVMapping(MultipleDVMapping mdv) {
		if (exd == null) throw new NullPointerException("External data set is NULL");
		
		DataSet ds = exd.getDataSet();
		if (ds == null) throw new NullPointerException("Dataset definition is NULL");
		
		Piecewise pw = mdv.getPiecewise();
		if (pw == null) return;
		
		Accessor a = lexer.getAccessor();
		ColumnDefinition col = a.fetchColumnDefinition(ds, mdv.getColumnRef());
		if (col == null) throw new NullPointerException("An MDV column reference is NULL.");
		
		List<Piece> pieces = pw.getListOfPiece();
		if (pieces.isEmpty()) return;
		
		for (int i = 0; i < pieces.size(); i++) {
			Piece piece = pieces.get(i);
			if (piece == null) continue;
			mdv_refs.add(new MultipleDvRef(lexer, col, piece));
		}
	}
}

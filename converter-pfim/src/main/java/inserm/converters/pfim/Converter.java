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

package inserm.converters.pfim;

import static crx.converter.engine.PharmMLTypeChecker.isColumnDefinition;
import static crx.converter.engine.PharmMLTypeChecker.isColumnMapping;
import static crx.converter.engine.PharmMLTypeChecker.isCommonParameter;
import static crx.converter.engine.PharmMLTypeChecker.isContinuousCovariate;
import static crx.converter.engine.PharmMLTypeChecker.isCovariate;
import static crx.converter.engine.PharmMLTypeChecker.isDerivative;
import static crx.converter.engine.PharmMLTypeChecker.isFunction;
import static crx.converter.engine.PharmMLTypeChecker.isFunctionCall;
import static crx.converter.engine.PharmMLTypeChecker.isFunctionParameter;
import static crx.converter.engine.PharmMLTypeChecker.isGeneralError;
import static crx.converter.engine.PharmMLTypeChecker.isIndependentVariable;
import static crx.converter.engine.PharmMLTypeChecker.isIndividualParameter;
import static crx.converter.engine.PharmMLTypeChecker.isInt;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isLogicalBinaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isObservationError;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;
import static crx.converter.engine.PharmMLTypeChecker.isRandomVariable;
import static crx.converter.engine.PharmMLTypeChecker.isStructuredError;
import static crx.converter.engine.PharmMLTypeChecker.isSymbol;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isVariabilityLevelDefinition;
import static crx.converter.engine.PharmMLTypeChecker.isVariableReference;
import static crx.converter.engine.assoc.DependencyRef.createElementsUnderConsideration;
import static crx.converter.engine.assoc.DependencyRef.updateDependencyContext;
import static eu.ddmore.libpharmml.dom.dataset.ColumnType.ADM;
import static eu.ddmore.libpharmml.dom.dataset.ColumnType.DOSE;
import inserm.converters.pfim.parts.CovariateBlockImpl;
import inserm.converters.pfim.parts.ObservationBlockImpl;
import inserm.converters.pfim.parts.OptimalDesignStepImpl;
import inserm.converters.pfim.parts.ParameterBlockImpl;
import inserm.converters.pfim.parts.StructuralBlockImpl;
import inserm.converters.pfim.parts.TrialDesignBlockImpl;
import inserm.converters.pfim.tree.TreeMaker_;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.Accessor;
import crx.converter.engine.ConversionDetail_;
import crx.converter.engine.Manager;
import crx.converter.engine.ParameterContext;
import crx.converter.engine.Part;
import crx.converter.engine.ScriptDefinition;
import crx.converter.engine.SymbolReader;
import crx.converter.engine.Utils;
import crx.converter.engine.VariableDeclarationContext;
import crx.converter.engine.assoc.DependencyGraph;
import crx.converter.engine.assoc.DependencyLexer;
import crx.converter.engine.assoc.DependencyRef;
import crx.converter.engine.common.DataFiles;
import crx.converter.engine.common.ObservationParameter;
import crx.converter.engine.common.SimulationOutput;
import crx.converter.engine.common.TabularDataset;
import crx.converter.spi.IParser;
import crx.converter.spi.OptimalDesignLexer;
import crx.converter.spi.blocks.CovariateBlock;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.blocks.TrialDesignBlock;
import crx.converter.spi.blocks.TrialDesignBlock2;
import crx.converter.spi.blocks.VariabilityBlock;
import crx.converter.spi.steps.EstimationStep;
import crx.converter.spi.steps.OptimalDesignStep_;
import crx.converter.spi.steps.SimulationStep;
import crx.converter.tree.BaseTreeMaker;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.NestedTreeRef;
import crx.converter.tree.Node;
import crx.converter.tree.TreeMaker;
import eu.ddmore.convertertoolbox.api.domain.LanguageVersion;
import eu.ddmore.convertertoolbox.api.domain.Version;
import eu.ddmore.convertertoolbox.api.response.ConversionDetail;
import eu.ddmore.convertertoolbox.api.response.ConversionReport;
import eu.ddmore.convertertoolbox.api.response.ConversionReport.ConversionCode;
import eu.ddmore.convertertoolbox.domain.ConversionReportImpl;
import eu.ddmore.convertertoolbox.domain.LanguageVersionImpl;
import eu.ddmore.convertertoolbox.domain.VersionImpl;
import eu.ddmore.libpharmml.ILibPharmML;
import eu.ddmore.libpharmml.IPharmMLResource;
import eu.ddmore.libpharmml.IValidationReport;
import eu.ddmore.libpharmml.PharmMlFactory;
import eu.ddmore.libpharmml.dom.IndependentVariable;
import eu.ddmore.libpharmml.dom.PharmML;
import eu.ddmore.libpharmml.dom.commontypes.CommonVariableDefinition;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.FunctionDefinition;
import eu.ddmore.libpharmml.dom.commontypes.FunctionParameter;
import eu.ddmore.libpharmml.dom.commontypes.IntValue;
import eu.ddmore.libpharmml.dom.commontypes.Name;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLElement;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.Symbol;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.dataset.ColumnDefinition;
import eu.ddmore.libpharmml.dom.dataset.ColumnMapping;
import eu.ddmore.libpharmml.dom.dataset.ColumnReference;
import eu.ddmore.libpharmml.dom.dataset.ColumnType;
import eu.ddmore.libpharmml.dom.dataset.DataSet;
import eu.ddmore.libpharmml.dom.dataset.HeaderColumnsDefinition;
import eu.ddmore.libpharmml.dom.dataset.MapType;
import eu.ddmore.libpharmml.dom.dataset.TargetMapping;
import eu.ddmore.libpharmml.dom.maths.Condition;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType.FunctionArgument;
import eu.ddmore.libpharmml.dom.maths.LogicBinOp;
import eu.ddmore.libpharmml.dom.maths.Piece;
import eu.ddmore.libpharmml.dom.maths.Piecewise;
import eu.ddmore.libpharmml.dom.modeldefn.CommonParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ContinuousCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.ContinuousObservationModel;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateModel;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateRelation;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateTransformation;
import eu.ddmore.libpharmml.dom.modeldefn.FixedEffectRelation;
import eu.ddmore.libpharmml.dom.modeldefn.GeneralObsError;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ModelDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.ObservationError;
import eu.ddmore.libpharmml.dom.modeldefn.ObservationModel;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterModel;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomEffect;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modeldefn.StructuralModel;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel.LinearCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError;
import eu.ddmore.libpharmml.dom.modeldefn.TransformedCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.VariabilityDefnBlock;
import eu.ddmore.libpharmml.dom.modeldefn.VariabilityLevelDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.pkmacro.PKMacro;
import eu.ddmore.libpharmml.dom.modeldefn.pkmacro.PKMacroList;
import eu.ddmore.libpharmml.dom.modellingsteps.ModellingSteps;
import eu.ddmore.libpharmml.dom.modellingsteps.OptimalDesignStep;
import eu.ddmore.libpharmml.dom.trialdesign.ExternalDataSet;
import eu.ddmore.libpharmml.dom.trialdesign.TrialDesign;
import eu.ddmore.libpharmml.dom.uncertml.VarRefType;
import eu.ddmore.libpharmml.impl.PharmMLVersion;
import eu.ddmore.libpharmml.pkmacro.exceptions.InvalidMacroException;
import eu.ddmore.libpharmml.pkmacro.translation.Input;
import eu.ddmore.libpharmml.pkmacro.translation.MacroOutput;
import eu.ddmore.libpharmml.pkmacro.translation.Translator;

/**
 * PFIM converter.
 */
public class Converter extends DependencyLexer implements OptimalDesignLexer {
	private static eu.ddmore.libpharmml.dom.commontypes.ObjectFactory cof = new eu.ddmore.libpharmml.dom.commontypes.ObjectFactory();
	private static eu.ddmore.libpharmml.dom.dataset.ObjectFactory ds_of = new eu.ddmore.libpharmml.dom.dataset.ObjectFactory();
	private static Manager m = new Manager();	
	private static final String MAJOR_VERSION = "4";
	private static eu.ddmore.libpharmml.dom.maths.ObjectFactory math_of = new eu.ddmore.libpharmml.dom.maths.ObjectFactory();
	/**
	 * Converter name
	 */
	public static final String NAME = "PFIM";
	
	private static Piece addPiece(Piecewise block, boolean createCondition) {
		Piece piece = null;
		
		if (block != null) {
			piece = new Piece();
			block.getListOfPiece().add(piece);
			
			if (createCondition) {
				Condition cond = new Condition();
				piece.setCondition(cond);
			}
		}
		
		return piece;
	}
	
	private static JAXBElement<IntValue> createScalar(Integer value) {
		IntValue o = new IntValue(value);
		return cof.createInt(o);
	}
	
	public static LogicBinOp logicalBinaryOperator(String pharmml_symbol, Object [] elements, Accessor a) {
		if (elements == null || pharmml_symbol == null || a == null) 
			throw new NullPointerException("A function argument cannot be NULL."); 
		
		LogicBinOp op = new LogicBinOp();
		op.setOp(pharmml_symbol);
	
		int addedElements = 0;
		for (Object o : elements) {
			if (isInteger(o)) {
				op.getContent().add(createScalar((Integer) o));
				addedElements++;
			} else if (isLogicalBinaryOperation(o)) {
				op.getContent().add(math_of.createLogicBinop((LogicBinOp) o));
				addedElements++;
			} else if (isColumnDefinition(o)) {
				ColumnDefinition col = (ColumnDefinition) o;
				ColumnReference ref = new ColumnReference();	
				ref.setColumnIdRef(col.getColumnId());
				op.getContent().add(ds_of.createColumnRef(ref));
				addedElements++;
			}
			
		}
		if (addedElements != 2) throw new IllegalStateException("A binary operator can only contains 2 elements."); 
		
		return op;
	}
	
	private static boolean replace(PharmML dom, StructuralModel old_sm, StructuralModel new_sm) {
		if (dom == null || old_sm == null || new_sm == null) return false;
		
		String blkId = old_sm.getBlkId();
		if (blkId == null) return false;
		
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return false;
		
		List<StructuralModel> sms_original = def.getListOfStructuralModel();
		if (sms_original == null) return false;
		if (sms_original.isEmpty()) return false;
		
		List<StructuralModel> sms = new ArrayList<StructuralModel>();
		sms.addAll(sms_original);
		
		boolean foundBlock = false;
		for (StructuralModel sm : sms) {
			if (sm == null) continue;
			if (sm.equals(old_sm) ) {
				foundBlock = true;
				break;
			}
		}
		
		if (!foundBlock) return false;
		
		// Replace the old structural model with the new one.
		new_sm.setBlkId(blkId);
		sms_original.clear();
		for (StructuralModel sm : sms) {
			if (sm == null) continue;
			if (sm.equals(old_sm) ) sms_original.add(new_sm);
			else sms_original.add(sm);
		}
		
		return true;
	}
	
	public static SymbolRef symbolRef(PharmMLRootType o, Accessor a) {
		String symbId = null;
		
		boolean addScope = false;
		
		if (isSymbolReference(o)) {
			return (SymbolRef) o;
		} else if (isCommonParameter(o)) { 
			symbId = ((CommonParameter) o).getSymbId();
			addScope = true;
		} else if (isLocalVariable(o)) {
			symbId = ((VariableDefinition) o).getSymbId();
			addScope = true;
		} else if (isDerivative(o)) { 
			symbId = ((DerivativeVariable) o).getSymbId();
			addScope = true;
		} else if (isIndividualParameter(o)) { 
			symbId = ((IndividualParameter) o).getSymbId();
			addScope = true;
		} else if (isRandomVariable(o)) {
			symbId = ((ParameterRandomVariable) o).getSymbId();
			addScope = true;
		} else if (isIndependentVariable(o)) {
			symbId = ((IndependentVariable) o).getSymbId();
		} else if (isCovariate(o)) {
			symbId = ((CovariateDefinition) o).getSymbId();
			addScope = true;
		} else if (isFunctionParameter(o)) {
			symbId = ((FunctionParameter) o).getSymbId();
		} else if (isFunction(o)) {
			symbId = ((FunctionDefinition) o).getSymbId();
		} else if (isObservationError(o)) {
			symbId = ((ObservationError) o).getSymbId();
			addScope = true;
		} else if (isColumnDefinition(o)) {
			symbId = ((ColumnDefinition) o).getColumnId();
			addScope = false;	
		} else if (isContinuousCovariate(o)) {
			ContinuousCovariate ccov = (ContinuousCovariate) o;
			
			// INFO: Assuming a unary application for this release. 
			for (CovariateTransformation trans : ccov.getListOfTransformation()) {
				if (trans == null) continue;
				
				TransformedCovariate tc = trans.getTransformedCovariate();
				if (tc == null) continue;
				
				symbId = tc.getSymbId();
				addScope = true;
				break;
			}
		} else if (isVariabilityLevelDefinition(o)) {
			VariabilityLevelDefinition level = (VariabilityLevelDefinition) o;
			symbId = level.getSymbId();
			addScope = true;
		}
		else if (isGeneralError(o)) {
			GeneralObsError goe = (GeneralObsError) o;
			symbId = goe.getSymbId();
			addScope = true;
		}
		else 
			throw new UnsupportedOperationException("Unsupported Symbol reference (src='" + o + "')");
		
		if (symbId == null) throw new NullPointerException("SymbId is NULL.");
		
		SymbolRef ref = new SymbolRef();
		ref.setSymbIdRef(symbId);
		
		if (addScope) {
			String blkId = a.getBlockId(o);
			if (blkId == null) {
				throw new NullPointerException("BlkId is not known (symbId='" + symbId + "', class='" + Utils.getClassName(o) + "')");
			}
			
			ref.setBlkIdRef(blkId);
		}
		
		return ref;
	}
	
	private boolean add_plotting_block = false;	
	private List<ConversionDetail_> cached_details = new ArrayList<ConversionDetail_>();
	private Version converterVersion = new VersionImpl(1, 7, 0);
	private boolean created_parameter_context = false;
	private StructuralBlock currentSb = null;
	private DataFiles data_files = new DataFiles();
	private boolean exceptionWithContinuousCovariate = false;
	private boolean hasResetColumnUsageRegToCov = false;
	private boolean is_echo_exception = true;
	private boolean isolate_dt = true;
	private boolean lexOnly = false;
	private ILibPharmML lib = null;
	private Map<StructuralModel, List<PKMacro>> macro_input_map = new HashMap<StructuralModel, List<PKMacro>>();
	private Map<StructuralModel, MacroOutput> macro_output_map = new HashMap<StructuralModel, MacroOutput>();
	private String model_filename = new String();
	private String name = "Core Lexer";
	private OptimalDesignStepImpl od_step = null;
	private List<String> ordered_levels = new ArrayList<String>();
	private String output_dir = ".";
	private Map<PopulationParameter, ParameterContext> param_context_map = new HashMap<PopulationParameter, ParameterContext>();
	private Parser parser;
	private ParameterBlock pb = null;
	private boolean permitEmptyTrialDesignBlock = false;
	private boolean remove_illegal_chars = false, filter_reserved_words = false;
	private IPharmMLResource res = null;
	private boolean resetRegToCov = true;
	private String run_id = null;
	private boolean save_renamed_symbol_list = true;
	private ScriptDefinition sd = new ScriptDefinition();
	private boolean sort_structural_model = true;
	
	/**
	 * Language source
	 */
	protected LanguageVersion source = null;
	
	/**
	 * Target language
	 */
	protected LanguageVersion target = null;
	private PharmMLVersion target_level = PharmMLVersion.V0_8_1;
	private Translator tr = new Translator();
	private boolean translate_macros = true;
	private boolean usePiecewiseAsEvents = false;
	private IValidationReport validation_report = null;
	
	public Converter() throws IOException, NullPointerException {
		super();
		initLibrary();
		
		setUsePiecewiseAsEvents(true);
		tm = new TreeMaker_();
		name = NAME;
		
		// Register lexer/parser dependency.
		Parser p = new Parser();
		setParser(p);
		p.setLexer(this);
		
		source = createPharmMLVersion();
		
		VersionImpl target_version = new VersionImpl(Integer.parseInt(MAJOR_VERSION), 0, 0);
		target = new LanguageVersionImpl("PFIM", target_version);
	}
	
	@Override
	public void addIndexSymbol(Object key, String value) { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean addStatement(NestedTreeRef ref) {
		if (ref != null) {
			if (ref.element != null && ref.bt != null) {
				if (sd.getStatementsMap().containsKey(ref.element)) sd.getStatementsMap().remove(ref.element);
				sd.getStatementsMap().put(ref.element, ref.bt);
				return true;
			}
		}
		
		return false;
	}
	
	@Override
	public void addStatement(Object element, BinaryTree tree) {
		if (element != null && tree != null) sd.getStatementsMap().put(element, tree);
	}
	
	private void buildPartTrees(Map<String, Part> parts) {
		if (parts == null) return;
		Set<String> keys = parts.keySet();
		for (String key : keys) {
			if (key == null) continue;
			Part part = parts.get(key);
			if (part == null) continue;
			part.buildTrees();
		}
	}
	
	private void checkForContinuousCovariates() {
		for (CovariateBlock cb : getCovariateBlocks()) {
			if (cb.getContinuousCovariates().size() > 0) {
				if (exceptionWithContinuousCovariate) throw new IllegalStateException("Input PharmML model continuous covariate. These are not supported by PFIM.");
				else {
					System.err.println("WARNING: Input model has continuous covariates. These are not supported by PFIM.");
					ConversionDetail_ detail = new ConversionDetail_();
					detail.addInfo("warning", "Input model has continuous covariates. These are not supported by PFIM.");
					detail.setSeverity(ConversionDetail.Severity.WARNING);
					cached_details.add(detail);
				}
			}
		}
	}
	
	@Override
	public void checkForTrialDesignIfRandomisedModel() {
		boolean has_mixed_effect = false;
		for (StructuralBlock sb : getStructuralBlocks()) {
			if (sb.isMixedEffect()) {
				has_mixed_effect = true;
				break;
			}
		}
		
		if (has_mixed_effect) {
			if (getTrialDesign() == null) 
				throw new IllegalStateException("Trial design not included in model even though contains individual parameter terms.");
		}
	}
	
	/**
	 * Parse the blocks in a PharmML model.
	 * @param outputDirectory Output Directory
	 * @throws IOException
	 */
	public void createBlocks(File outputDirectory) throws IOException { 
		parser.initialise();
		if (outputDirectory == null) throw new NullPointerException("The output directroy is NULL");
		
		createFunctionTrees();
		createVariabilityMap();
		createObservationMap();
		createTrialDesignMaps();
		createParameterMap();
		createStepMap();
		createStepTrees();
		createStateMap();
		createCovariateMap();
		
		createParameterTrees();
		createStructuralTrees();
		createCovariateTrees();
		
		checkForContinuousCovariates();
		sortElementOrdering();
		
		ParameterBlockImpl pb_ = (ParameterBlockImpl) getParameterBlock();
		if (pb_ != null) pb_.buildFixedEffectRefs();
	}
	
	private void createCovariateMap() {
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<CovariateModel> cmts = def.getListOfCovariateModel();
		if (cmts == null) return;
		
		for (CovariateModel cmt : cmts) getCovariateBlocks().add(new CovariateBlockImpl(cmt, this));
	}
	
	private void createCovariateTrees() { for (CovariateBlock cb : getCovariateBlocks()) cb.buildTrees(); }
	
	/**
	 * Pre-code generation parsing of any dose settings mapped to a input data column
	 * in an external data set.
	 */
	private void createExternalFileDoseVariableMappings() {
		if (!translate_macros) return;
		getAccessor();
		
		TrialDesign td = dom.getTrialDesign();
		if (td == null) return;
		
		List<ExternalDataSet> eds = td.getListOfExternalDataSet();
		if (eds == null) return;
		
		for (ExternalDataSet ed : eds) {
			if (ed == null) continue;
			
			List<PharmMLRootType> col_elements = ed.getListOfColumnMappingOrColumnTransformationOrMultipleDVMapping();
			DataSet ds = ed.getDataSet();
			
			if (col_elements == null || ds == null) continue;
			
			for (PharmMLRootType col_element : col_elements) {
				if (isColumnMapping(col_element)) {
					ColumnMapping mapping = (ColumnMapping) col_element;
					ColumnDefinition col = accessor.fetchColumnDefinition(ds, mapping.getColumnRef());
					if (col == null) continue;
					
					ColumnType usage = col.getListOfColumnType().get(0);
					if (usage == null) continue;
					else if (usage.equals(ColumnType.DOSE)) {
						if (mapping.getPiecewise() != null) continue;
					}
				}
			}
		}
	}
	
	/**
	 * Create the trees representing model functions.
	 */
	private void createFunctionTrees() {
		if (dom == null) return;
		List<FunctionDefinition> funcs = dom.getListOfFunctionDefinition();
		
		if (funcs == null) return;
		else if (funcs.isEmpty()) return;

		for (FunctionDefinition func : funcs) {
			BinaryTree bt = tm.newInstance(func);
			if (bt == null) throw new IllegalStateException("The funcion assignment block cannot be null.");
			sd.getStatementsMap().put(func, bt);
			sd.getFunctions().add(func);
		}
	}
	
	private void createObservationMap() {
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<ObservationModel> observations = def.getListOfObservationModel();
		for (ObservationModel ob : observations) {
			if (ob == null) continue;
			ObservationBlock obb = new ObservationBlockImpl(ob, this);
			obb.buildTrees();
			getObservationBlocks().add(obb);
		}
	}
	
	/**
	 * Create a parameter block
	 * @param pmt Parent parameter block
	 * @return ParameterBlock
	 */
	private ParameterBlock createParameterBlock(ParameterModel pmt) { return new ParameterBlockImpl(pmt, this); }
	
	@Override
	public void createParameterContext() {
		if (created_parameter_context) return;

		ParameterBlock pdb = getParameterBlock();
		if (pdb == null) return;

		List<PopulationParameter> params = pdb.getParameters();
		for (PopulationParameter p : params) getParameterContext(p);

		doParameterContext_ObservationBlocks();
		doParameterContext_Individuals(pdb);

		created_parameter_context = true; // So not to call this method twice.
	}
	
	private void createParameterMap() {
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<ParameterModel> param_models = def.getListOfParameterModel();
		if (param_models.isEmpty()) return;
		
		for (ParameterModel pm : param_models)
			if (pm == null) continue;
			else sd.getParameterBlocks().add(new ParameterBlockImpl(pm, this));
	}
	
	private void createParameterTrees() {
		for (ParameterBlock pb : sd.getParameterBlocks()) {
			if (pb == null) continue;
			pb.buildTrees();
		}
	}
	
	/**
	 * Create the language version/level handle for the preferred PharmML release.
	 * @return LanguageVersion
	 */
	protected LanguageVersion createPharmMLVersion() {
		VersionImpl preferred_pharmml_version = new VersionImpl(0, 8, 1);
		return new LanguageVersionImpl("PharmML", preferred_pharmml_version);
	}
	
	protected File createScript(File src, File outputDirectory) throws Exception {
        if (outputDirectory == null) throw new NullPointerException("The output directroy is NULL");
        String output_filename = parser.getScriptFilename(outputDirectory.getAbsolutePath());
        
        parser.writeSTDIN();
        parser.writeBatchFile();
        
        PrintWriter fout = new PrintWriter(output_filename);
	
        parser.writePreMainBlockElements(fout, outputDirectory);
        parser.writeModelCall(fout);
        parser.writeEOF(fout);
        fout.close();
        fout = null;
        parser.cleanUp();
        
        if (save_renamed_symbol_list) {
			SymbolReader z = parser.getSymbolReader();
			if (z != null) {
				if (z.hasModifiedSymbols()) z.saveChanges(getOutputDirectory());
			}
		}
        
        File f = new File(output_filename);
        return f; 
    }
	
	private void createStateMap() {
		if (dom == null) return;
		
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<StructuralModel> smts = def.getListOfStructuralModel();
		if (smts.isEmpty()) return;
	
		for (StructuralModel sm : smts) 
			if (sm  == null) continue;
			else {
				StructuralBlock sb = createStructuralBlock(sm);
				if (macro_input_map.containsKey(sm)) sb.setPKMacros(macro_input_map.get(sm));
				if (macro_output_map.containsKey(sm)) sb.setPKMacroOutput(macro_output_map.get(sm)); 
				getStructuralBlocks().add(sb);
			}
	}
	
	/**
	 * Create the step map of registered tasks to be done to a model.
	 */
	private void createStepMap() {
		if (dom == null) return;
		
		ModellingSteps steps_ = dom.getModellingSteps();
		TrialDesign td = dom.getTrialDesign();
		
		if (steps_ == null) throw new IllegalStateException("The PharmML file does not have the required 'ModellingSteps' section.");
		if (td == null) throw new NullPointerException("Model trial design is not specified.");
		
		List<OptimalDesignStep> od_steps = steps_.getListOfOptimalDesignStep();
		if (od_steps == null) throw new NullPointerException("Optimal design task not found");
		if (od_steps.isEmpty()) throw new IllegalStateException("An optimal design step not specified");
	
		od_step = new OptimalDesignStepImpl(od_steps.get(0), this);
		sd.getStepsMap().put(od_step.toString(), od_step);
	}
	
	private void createStepTrees() {  buildPartTrees(sd.getStepsMap()); }
	
	private StructuralBlock createStructuralBlock(StructuralModel sm) { return new StructuralBlockImpl(sm, this); }
	
	private void createStructuralTrees() {
		for (StructuralBlock sb : getStructuralBlocks()) {
			if (sb == null) continue;
			sb.buildTrees();
		}
	}
	
	@Override
	public BinaryTree createTree(Object o) {
		BinaryTree bt = tm.newInstance(o);
		addStatement(o, bt);
		updateNestedTrees();
		return bt;
	}
	
	/**
	 * Factory method to create a trial design block.
	 * The trial design block does basic processing of a trial design schema in PharmML.
	 * If the default processing of TD is not appropriate, override this method and replace
	 * with different implementation of the TrialDesign block.
	 * 
	 * @param td Trial Design Model
	 * @return TrialDesignBlock
	 */
	private TrialDesignBlock createTrialDesignBlock(TrialDesign td) { return new TrialDesignBlockImpl(td, this); }
	
	private void createTrialDesignMaps() throws IOException {
		if (dom == null) return;
		
		TrialDesign td = dom.getTrialDesign();
		
		if (td == null) return; 
		else {
			TrialDesignBlock tdb = createTrialDesignBlock(td);
			tdb.buildTrees();
			
			sd.setTrialDesignBlock(tdb);
		}
		
		sd.getBlocksMap().put(getTrialDesign().getName(), getTrialDesign()); // Use the same event structure as step classes.
	}
	
	private void createVariabilityMap() {
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<VariabilityDefnBlock> vms = def.getListOfVariabilityModel();
		if (vms == null) return;
		
		for (VariabilityDefnBlock vm : vms) {
			if (vm == null) continue;
			VariabilityBlock vb = new inserm.converters.pfim.parts.VariabilityBlockImpl(vm, this);
			sd.getVariabilityBlocks().add(vb);
		}
	}
	
	private void doParameterContext_CovariateRelation(List<CovariateRelation> covs) {
		if (covs == null) return;
		if (covs.isEmpty()) return;

		for (CovariateRelation cov : covs) {
			if (cov == null) continue;
			for (FixedEffectRelation fixed_effect : cov.getListOfFixedEffect()) {
				if (fixed_effect == null) continue;
				SymbolRef ref = fixed_effect.getSymbRef();
				if (ref == null) continue;
				Object o = accessor.fetchElement(ref);
				if (isPopulationParameter(o)) {
					PopulationParameter p = (PopulationParameter) o;
					ParameterContext ctx = getParameterContext(p);
					if (ctx != null) ctx.theta_fixed_effect = true;
				}
			}
		}
	}
	
	private void doParameterContext_Individuals(ParameterBlock pdb) {
		// Look for 'Group' variable context, i.e. a generalised covariate.
		// Or a THETA context.
		List<IndividualParameter> ips = pdb.getIndividualParameters();
		for (IndividualParameter ip : ips) {
			if (ip == null) continue;

			if (ip.getStructuredModel() != null) {
				StructuredModel gm = ip.getStructuredModel();
				doParameterContext_LinearCovariate(gm.getLinearCovariate());
				doParameterContext_RandomEffects(gm.getListOfRandomEffects());
			}
		}
	}
	
	private void doParameterContext_LinearCovariate(LinearCovariate lc) {
		if (lc == null) return;
		
		doParameterContext_PopulationValue_(lc.getPopulationValue());
		doParameterContext_CovariateRelation(lc.getListOfCovariate());
	}

	private void doParameterContext_ObservationBlocks() {
		List<ObservationBlock> obs = getObservationBlocks();
		for (ObservationBlock ob : obs) {
			if (ob == null) continue;
			
			ObservationError error_model = ob.getObservationError();
			List<SymbolRef> residuals = new ArrayList<SymbolRef>();  
			List<NestedTreeRef> error_trees = new ArrayList<NestedTreeRef>();
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
				else
					residuals.add(residual_error.getSymbRef());
			} else if (isGeneralError(error_model)) {
				GeneralObsError goe = (GeneralObsError) error_model;
				error_trees.add(new NestedTreeRef(goe, tm.newInstance(goe)));
				error_trees.add(new NestedTreeRef(goe.getAssign(), tm.newInstance(goe.getAssign())));
			}
			
			for (NestedTreeRef ntr : error_trees) {
				if (ntr == null) continue;
				for (Node node : ntr.bt.nodes) {
					if (isFunctionCall(node.data)) {
						FunctionCallType call = (FunctionCallType) node.data;
						for (FunctionArgument arg : call.getListOfFunctionArgument()) {
							if (arg == null) continue;
							
							if (arg.getSymbRef() == null) continue; // Native type so 2 point squinting at the type.
							
							PharmMLRootType element = accessor.fetchElement(arg.getSymbRef());
							if (isPopulationParameter(element)) {
								ParameterContext ctx = getParameterContext((PopulationParameter) element);
								if (ctx != null) ctx.error_model = true;
							}
						}
					} 
				}
			}
		}
	}

	private void doParameterContext_PopulationValue_(StructuredModel.LinearCovariate.PopulationValue pv) {
		if (pv == null) return;
		BinaryTree bt = tm.newInstance(pv);
		
		List<BinaryTree> bts = new ArrayList<BinaryTree>();
		bts.add(bt);
		for (NestedTreeRef ref : tm.getNestedTrees()) {
			if (ref == null) continue;
			bts.add(bt);
		}
		
		for (BinaryTree bt_ : bts) {
			for (Node node : bt_.nodes) {
				if (node == null) continue;
				if (isSymbolReference(node.data)) {
					SymbolRef sref = (SymbolRef) node.data;
					PharmMLRootType element = accessor.fetchElement(sref);
					if (isPopulationParameter(element)) {
						PopulationParameter p = (PopulationParameter) element;
						ParameterContext ctx = getParameterContext(p);
						if (ctx == null) continue;
						ctx.theta = true;
					}
				}
			}
		}
	}
	
	private void doParameterContext_RandomEffects(List<ParameterRandomEffect> rfs) {
		if (rfs == null) return;
		if (rfs.isEmpty()) return;

		for (ParameterRandomEffect rf : rfs) {
			if (rf == null) continue;
			for (SymbolRef ref : rf.getSymbRef()) {
				if (ref == null) continue;
				PharmMLRootType element = accessor.fetchElement(ref);
				if (element != null) {
					List<BinaryTree> bts = new ArrayList<BinaryTree>();
					bts.add(tm.newInstance(element));
					for (NestedTreeRef nref : tm.getNestedTrees()) {
						if (ref == null) continue;
						bts.add(nref.bt);
					}

					for (BinaryTree bt : bts) {
						for (Node node : bt.nodes) {
							if (isVariableReference(node.data)) {
								VarRefType vref = (VarRefType) node.data; 
								PharmMLRootType o = accessor.fetchElement(vref);
								if (isPopulationParameter(o)) {
									PopulationParameter p = (PopulationParameter) o;
									ParameterContext ctx = getParameterContext(p);
									if (ctx == null) continue;
									ctx.omega = true;
								}
							} else if (isSymbolReference(node.data)) {
								SymbolRef sref = (SymbolRef) node.data;
								PharmMLRootType o = accessor.fetchElement(sref);
								if (isPopulationParameter(o)) {
									PopulationParameter p = (PopulationParameter) o;
									ParameterContext ctx = getParameterContext(p);
									if (ctx == null) continue;
									ctx.omega = true;
								}
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Access the accessor handle
	 * @return Accessor
	 */
	public Accessor getAccessor() {
		if (accessor == null) accessor = new Accessor(dom); 
		return accessor;
	}
	
	@Override
	public Version getConverterVersion() { return converterVersion; }
	
	@Override
	public CovariateBlock getCovariateBlock() {
		List<CovariateBlock> cbs = getCovariateBlocks();
		if (cbs != null) {
			if (!cbs.isEmpty()) return cbs.get(0);
		}
		
		return null;
	}
	
	@Override
	public List<CovariateBlock> getCovariateBlocks() {
		return sd.getCovariateBlocks();
	}
	
	@Override
	public List<CovariateDefinition> getCovariates() {
		List<CovariateDefinition> list = new ArrayList<CovariateDefinition>();
		
		for (CovariateBlock cb : getCovariateBlocks()) {
			List<CovariateDefinition> covs = cb.getCovariates();
			for (CovariateDefinition cov : covs) {
				if (cov == null) continue;
				list.add(cov);
			}
		}
		
		return list;
	}
	
	/**
	 * Generate a conversion report.
	 * @param path Path to the generated file.
	 * @return eu.ddmore.convertertoolbox.api.response.ConversionReport
	 */
	private ConversionReport getCrxSuccessReport(File path) {
		ConversionReport report = new ConversionReportImpl();
		ConversionDetail_ detail = new ConversionDetail_();
		detail.addInfo("status", "script created.");
		if (path != null) detail.setFile(path);
		detail.setSeverity(ConversionDetail.Severity.INFO);
		report.addDetail(detail);
		report.setReturnCode(ConversionCode.SUCCESS);
		
		for (ConversionDetail_ cached_detail : cached_details) report.addDetail(cached_detail);
		
		return report;
	}
	
	@Override
	public DataFiles getDataFiles() {
		return data_files;
	}
	
	@Override
	public List<ObservationBlock> getErrorModels() {
		return getObservationBlocks();
	}
	
	@Override
	public EstimationStep getEstimationStep() { throw new UnsupportedOperationException(); }
	
	/**
	 * Generate a conversion report with an exception
	 * @param e Exception
	 * @return eu.ddmore.convertertoolbox.api.response.ConversionReport
	 */
	private ConversionReport getExceptionReport(Exception e) {
		ConversionReport report = new ConversionReportImpl();
		
		if (is_echo_exception) e.printStackTrace(System.err);
		
		ConversionDetail_ detail = new ConversionDetail_();
		detail.addInfo("error_message", e.getMessage());
		detail.setSeverity(ConversionDetail.Severity.ERROR);
		detail.setMessage(e.getMessage());
		report.addDetail(detail);
		report.setReturnCode(ConversionCode.FAILURE);
		
		for (ConversionDetail_ cached_detail : cached_details) report.addDetail(cached_detail);
		
		return report;
	}
	
	/**
	 * Get list of exported variables associated with observation models.
	 * @return List<VariableDefinition>
	 */
	public List<VariableDefinition> getExportedLocalVariables() {
		List<VariableDefinition> exported_variables = new ArrayList<VariableDefinition>();
		for (ObservationBlock ob : getObservationBlocks()) {
			if (ob == null) continue;
			
			for (SimulationOutput output : ob.getSimulationOutputs()) 
				if (isLocalVariable(output.v)) exported_variables.add((VariableDefinition) output.v);
		}
		
		return exported_variables;
	}
	
	
	
	@Override
	public String getIndexSymbol(Object key) { throw new UnsupportedOperationException(); }
	
	@Override
	public Integer getIndividualParameterIndex(String symbol) { throw new UnsupportedOperationException(); }
	
	@Override
	// Return zero index by default.
	public List<VariableDefinition> getLocalVariables() {
		if (!sd.getStructuralBlocks().isEmpty()) {
			StructuralBlock structuralBlock= getStructuralBlocks().get(0);
			if (structuralBlock != null) return structuralBlock.getLocalVariables();
			else return new ArrayList<VariableDefinition>();
		} else 
			return new ArrayList<VariableDefinition>();
	}
	
	@Override
	public String getModelFilename() { return model_filename; }
	
	@Override
	public String getModelName() {
		String name = null;
		
		if (dom != null) {
			Name n = dom.getName();
			if (n == null) return name;
			else {
				name = n.getValue();
				if (name != null) {
					name = name.replaceAll("\n", "");
					name = name.replaceAll("\\s+", " ");
					name = name.replace('\\', '/'); 
				}
			}		
		}
		
		return name;
	}
	
	@Override
	public Integer getModelParameterIndex(String name) {
		Integer idx = -1;
		if (!sd.getParameterBlocks().isEmpty()) {
			if (getParameterBlock() != null) idx = getParameterBlock().getParameterIndex(name);
		}
		return idx;
	}
	
	@Override
	// Return zero index by default.
	public List<PopulationParameter> getModelParameters() {
		if (!sd.getParameterBlocks().isEmpty()) {
			if (getParameterBlock() != null) return getParameterBlock().getParameters();
			else return new ArrayList<PopulationParameter>();
		} else 
			return new ArrayList<PopulationParameter>();
	}

	@Override
	public String getName() { return name; }
	
	@Override
	public List<ObservationBlock> getObservationBlocks() { return sd.getObservationBlocks(); }
	
	@Override
	public ObservationParameter getObservationParameter(PopulationParameter p) {
		ObservationParameter op = null;
		for (ObservationBlock ob : getObservationBlocks()) {
			if (ob == null) continue;
			if (ob.isObservationParameter(p)) {
				op = ob.getObservationParameter(p);
				break;
			}
		}
		
		return op;
	}
	
	@Override
	public OptimalDesignStep_ getOptimalDesignStep() { return od_step; }
	
	public OptimalDesignStepImpl getOptimalDesignStepRef() { return od_step; }
	
	@Override
	//  Overriding this method in case Java generates Windows style file paths that 'R' does not like.
	public String getOutputDirectory() {
		if (output_dir == null) return null;
		
		try {
			File dir = new File(output_dir);
			output_dir = dir.getCanonicalPath();
			output_dir = output_dir.replace('\\', '/');
		} catch (Exception e) {
			if (is_echo_exception) e.printStackTrace();
		}
		
		return output_dir; 
	}
	
	@Override
	public ParameterBlock getParameterBlock() {
		if (pb == null) {
			if (!sd.getParameterBlocks().isEmpty()) {
				pb = sd.getParameterBlocks().get(0); // Assume only 1 PM per model.
			} else {
				ParameterModel pmt = new ParameterModel();
				pmt.setBlkId("pm");
				ModelDefinition def = dom.getModelDefinition(); 
				if (def == null) {
					def = new ModelDefinition();
					dom.setModelDefinition(def);
				}
				def.getListOfParameterModel().add(pmt);
				pb = createParameterBlock(pmt);
			}
		}
			
		return pb;
	}
	
	private ParameterContext getParameterContext(PopulationParameter p) {
		ParameterContext ctx = null;
		
		if (p == null) throw new NullPointerException("Parameter is NULL.");
		if (param_context_map.containsKey(p)) ctx = param_context_map.get(p);
		else {
			ctx = new ParameterContext(p);
			param_context_map.put(p,  ctx);
		}
		
		return ctx;
	}
	
	@Override
	public Map<PopulationParameter, ParameterContext> getParameterContextMap() { return param_context_map; }

	@Override
	public IParser getParser() { return parser; }

	public ScriptDefinition getScriptDefinition() { return sd; }

	@Override
	public List<String> getSimulationOutputNames() {
		List<String> names = new ArrayList<String>();
		SimulationStep step = getSimulationStep();
		SymbolReader z = parser.getSymbolReader();
		
		if (step != null) {
			List<SymbolRef> outputs = step.getContinuousList();
			for (SymbolRef ref : outputs) {
				PharmMLRootType element = accessor.fetchElement(ref);
				if (element != null) names.add(z.get(element));
			}
		}
	
		return names;
	}

	
	
	@Override
	public Map<Integer, CommonVariableDefinition> getSimulationOutputs()  {
		StructuralBlock block = getStrucuturalBlock();
		return getSimulationOutputs(block);
	}
	
	@Override
	public Map<Integer, CommonVariableDefinition> getSimulationOutputs(StructuralBlock sb) { throw new UnsupportedOperationException(); }
	
	@Override
	public SimulationStep getSimulationStep() { throw new UnsupportedOperationException(); }
	
	@Override
	public List<String> getSortedVariabilityLevels() { return ordered_levels; }
	
	@Override
	public LanguageVersion getSource() { return source; }
	
	@Override
	public BinaryTree getStatement(Object key) {
		BinaryTree bt = null;
		
		if (key != null) {
			if (sd.getStatementsMap().containsKey(key)) {
				bt = sd.getStatementsMap().get(key);
				sd.getStatementsMap().remove(key);
			}
		}
		
		return bt;
	}
	
	@Override
	public Integer getStateVariableIndex(String name) {
		Integer idx = -1;
				
		PharmMLRootType element = accessor.fetchElement(name);
		if (element == null) return idx;
			
		StructuralBlock sb = getStrucuturalBlock();
		idx = sb.getStateVariableIndex(name);
					
		return idx;
	}
 
	@Override
	public List<StructuralBlock> getStructuralBlocks() { return sd.getStructuralBlocks(); }
	
	@Override
	public StructuralBlock getStrucuturalBlock() {
		if (currentSb == null) {
			if (!sd.getStructuralBlocks().isEmpty()) currentSb = sd.getStructuralBlocks().get(0);
		}
		
		return currentSb;
	}

	@Override
	public LanguageVersion getTarget() { return target; }
	
	private CommonVariableDefinition getTranslatedElement(Integer adm) {
		if (adm == null) return null;
		
		for (MacroOutput output : macro_output_map.values()) {
			if (output == null) continue;
			for (Input input : output.getListOfInput()) {
				if (input == null) continue;
				if (isInt(input.getAdm())) {
					IntValue i = (IntValue) input.getAdm();
					if (i.getValue() == null) continue;
					if (adm.intValue() == i.getValue().intValue()) return input.getTarget();
				}
			}
		}
		
		return null;
	}
	
	@Override
	public Translator getTranslator() { return tr; }
	
	@Override
	public TreeMaker getTreeMaker() {
		if (tm == null) {
			tm = new BaseTreeMaker();
			tm.setPermitDeclarationOnlyVariables(true);
		}
		
		return tm; 
	}
	
	@Override
	public TrialDesignBlock getTrialDesign() {
		return sd.getTrialDesignBlock();
	}
	
	@Override
	public IValidationReport getValidationReport() { return validation_report; }

	@Override
	public VariabilityBlock getVariabilityBlock(SymbolRef ref) {
		VariabilityBlock vb = null;
		
		if (ref != null) {
			String symbId = ref.getSymbIdRef();
			String blkId = ref.getBlkIdRef();
			
			if (blkId != null && symbId != null) {
				for (VariabilityBlock vb_ : sd.getVariabilityBlocks()) {
					if (vb_ == null) continue;
					if (blkId.equals(vb_.getName())) {
						if (vb_.hasSymbolId(symbId)) {
							vb = vb_;
							break;
						}
					}
				}
			}
			
			// Look-up with just the symbol name identifier if
			// the BlkId is not specified.
			if (vb == null) {
				if (symbId != null) {
					for (VariabilityBlock vb_ : sd.getVariabilityBlocks()) {
						if (vb_ == null) continue;
						if (vb_.hasSymbolId(symbId)) {
							vb = vb_;
							break;
						}
					}
				}
			}
		}
		
		return vb;
	}

	@Override
	public VariableDeclarationContext guessContext(VariableDefinition v) {
		if (v == null) return VariableDeclarationContext.UNKNOWN;
		
		if (v.getAssign() != null) return VariableDeclarationContext.ASSIGNED;
		
		String symbol = v.getSymbId();
		if (symbol == null) return VariableDeclarationContext.UNKNOWN;
		
		if (hasEstimation()) {
			EstimationStep est = getEstimationStep();
			if (est.hasTemporalDoseEvent()) {
				PharmMLRootType element = est.getTemporalDoseEvent().getTargetElement();
				if (isLocalVariable(element)) {
					VariableDefinition o = (VariableDefinition) element;
					if (symbol.equals(o.getSymbId())) {
						return VariableDeclarationContext.DT; 
					}
				}
			}
		}
		
		if (isTargetMappedDoseVariable(v)) return VariableDeclarationContext.DOSE;
		
		return VariableDeclarationContext.UNKNOWN;
	}
	
	@Override
	public boolean hasDoneEstimation() { throw new UnsupportedOperationException(); }

	@Override
	public boolean hasDosing() {
		if (getTrialDesign() != null) {
			TrialDesignBlockImpl o = (TrialDesignBlockImpl) getTrialDesign();
			return o.hasDosing();
		}
		
		return false;
	}
	
	@Override
	public boolean hasEstimation() {
		for (Part step : sd.getStepsMap().values()) if (step instanceof EstimationStep) return true;
		return false;
	}
	
	@Override
	public boolean hasEvents() {
		StructuralBlock sb = getStrucuturalBlock();
		if (sb != null) return sb.hasEvents();
		return false;
	}
	
	@Override
	public boolean hasExternalDatasets() { return !data_files.getExternalDataSets().isEmpty(); }
	
	@Override
	public boolean hasPlottingBlock() { return add_plotting_block; }

	/**
	 * Flag if the Lexer has adjustment an external dataset column usage declaration
	 * from a Monolix (Regressor) to a NONMEM (covariate) setting.
	 * @return boolean
	 */
	public boolean hasResetColumnFromRegToCov() { return hasResetColumnUsageRegToCov; }
	
	@Override
	public boolean hasSimulation() { throw new UnsupportedOperationException(); }

	@Override
	public boolean hasStatement(Object key) {
		if (key != null) return sd.getStatementsMap().containsKey(key);
		return false;
	}

	@Override
	public boolean hasTranslatedPKMacros() {
		if (!isTranslate()) return false;
		
		List<StructuralBlock> sbs = getStructuralBlocks();
		if (sbs == null) return false;
		
		for (StructuralBlock sb : sbs) {
			if (sb == null) continue;
			if (sb.getPKMacroOutput() != null) return true;
		}
		
		return false;
	}

	@Override 
	public boolean hasTrialDesign() { return getTrialDesign() != null; }

	@Override
	public boolean hasUntranslatedPKMacros() {
		if (isTranslate()) return false;
		
		List<StructuralBlock> sbs = getStructuralBlocks();
		if (sbs == null) return false;
		
		for (StructuralBlock sb : sbs) {
			if (sb == null) continue;
			if (sb.isUsingUntranslatedPKMacros()) return true;
		}
		
		return false;
	}
	
	@Override
	public boolean hasWashout() { throw new UnsupportedOperationException(); }
	
	private void initLibrary() { lib = PharmMlFactory.getInstance().createLibPharmML(); }
	
	@Override
	public boolean isADMEScript() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isAtLastStructuralBlock() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isCategoricalCovariate() {
		List<CovariateBlock> cbs = getCovariateBlocks();
		
		for (CovariateBlock cb : cbs) {
			if (cb == null) continue;
			if (cb.isCategorical()) return true;
		}
		
		return false;
	}
	
	@Override
	public boolean isCategoricalDiscrete() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isDDE() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isDiscrete() { return false; }
	
	@Override
	public boolean isDuplicateVariablesPermitted() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isFilterReservedWords() { return filter_reserved_words; }
	
	@Override
	public boolean isGeneratedDosingParameter(PopulationParameter p) { throw new UnsupportedOperationException();	 }
	
	@Override
	public boolean isIndexFromZero() { return true; }
	
	@Override
	public boolean isIndividualParameter_(String symbol) {
		if (symbol != null) {
			if (!sd.getParameterBlocks().isEmpty()) {
				if (getParameterBlock() != null) {
					List<IndividualParameter> ips = getParameterBlock().getIndividualParameters();
					for (IndividualParameter ip : ips) {
						if (ip == null) continue;
						String currentSymbol = ip.getSymbId();
						if (currentSymbol == null) continue;
						if (currentSymbol.equals(symbol)) return true;
					}
				}
			}
		}
		
		return false;
	}
	
	@Override
	public boolean isIsolateConditionalDoseVariable() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isIsolatingDoseTimingVariable() { return isolate_dt; }
	
	@Override
	public boolean isMixedEffect() { return false; }
	
	@Override
	public boolean isModelParameter(String name) {
		boolean isPopulationParameter = false;
		if (name == null) return isPopulationParameter; 
		if (!sd.getParameterBlocks().isEmpty()) {
			if (getParameterBlock() != null) isPopulationParameter = getParameterBlock().hasSymbolId(name);
		}
		
		return isPopulationParameter;
	}
	
	@Override
	public boolean isObservationParameter(PopulationParameter p) {
		if (p != null) {
			for (ObservationBlock ob : getObservationBlocks()) {
				if (ob == null) continue;
				if (ob.isObservationParameter(p)) return true;
			}
		}
		
		return false;
	}

	@Override
	public boolean isObservationParameter(SymbolRef ref) {
		if (ref != null) {
			Object element = accessor.fetchElement(ref);
			if (isPopulationParameter(element)) return isObservationParameter((PopulationParameter) element);
		}
		
		return false;
	}

	@Override
	public boolean isPermitEmptyTrialDesignBlock() { return permitEmptyTrialDesignBlock; }

	@Override
	public boolean isRemoveIllegalCharacters() { return remove_illegal_chars; }

	@Override
	public boolean isSaveSimulationOutput() { throw new UnsupportedOperationException();}
	
	@Override
	public boolean isStateVariable(String name) {
		boolean isState = false;
		
		if (!sd.getStructuralBlocks().isEmpty() && name != null) {
			for (StructuralBlock sb : sd.getStructuralBlocks()) {
				if (sb == null) continue;
				isState = sb.isStateVariable(name);
				if (isState) break;
			}
		}
		
		return isState;
	}
	
	@Override
	public boolean isStructuralBlockWithDosing(StructuralBlock sb) { throw new UnsupportedOperationException(); }
	
	private boolean isTargetMappedDoseVariable(VariableDefinition v) {
		if (v == null) return false;
		TrialDesignBlock2 tdb = (TrialDesignBlock2) getTrialDesign();
		if (tdb == null) return false;
		return tdb.getDoseTargetMap().containsValue(v);
	}
	
	@Override
	public boolean isTranslate() { return translate_macros; }
	
	@Override
	public boolean isTTE() {
		List<ObservationBlock> obs = getObservationBlocks();
		if (obs.isEmpty()) return false;
		
		for (ObservationBlock ob : obs) {
			if (ob == null) continue;
			if (ob.hasTimeToEventData()) return true;
		}
		
		return false;
	}
	
	@Override
	public boolean isUseGlobalConditionalDoseVariable() { throw new UnsupportedOperationException(); }
	
	@Override
	public boolean isUsePiecewiseAsEvents() { return usePiecewiseAsEvents; }
	
	@Override
	public boolean loadFunctionLibrary(File f)  { throw new UnsupportedOperationException(); }
	
	/**
	 * Load a PharmML file.
	 * @param xml_file_path File path to XML
	 * @throws IOException
	 */
	private void loadPharmML(File xml_file_path) throws IOException {
		InputStream in = new FileInputStream(xml_file_path);
		res = lib.createDomFromResource(in);
		in.close();
		
		dom = res.getDom();
		getAccessor();
		
		model_filename = xml_file_path.getAbsolutePath();
		parser.setPharmMLWrittenVersion(dom.getWrittenVersion());
		
		
		
		translatePKMacros();
		resetRegressorToCovariate();
		createExternalFileDoseVariableMappings();
	}
	
	@Override
    public ConversionReport performConvert(File src, File outputDirectory) {
        getScriptDefinition().flushAllSymbols();
        try {
        	if (run_id != null) parser.setRunId(run_id);
			else parser.setRunId(m.generateRunId());
            parser.setRunId(m.generateRunId());
            setOutputDirectory(outputDirectory);
            loadPharmML(src);
            
            // Parse each of the available structural blocks.
            createBlocks(outputDirectory);
            if (lexOnly) return getCrxSuccessReport(null);
            File f = createScript(src, outputDirectory);
            return getCrxSuccessReport(f);
        } catch (Exception e) {
            return getExceptionReport(e);
        }
    }
	
	@Override
	public void permitObjectiveETAs(boolean decision) { throw new UnsupportedOperationException(); }
	
	/**
	 * Re-map the dose target post PK Macro translation.
	 * Override as necessary 
	 * @param mapping Mapping Declaration
	 * @param table Column Looker-Upper
	 */
	private void remapColumnPostPKMacroTranslation(ColumnMapping mapping, TabularDataset table) {
		if (mapping == null) return;
		if (mapping.getListOfTargetMapping().isEmpty()) return;
		
		List<TargetMapping> targets = new ArrayList<TargetMapping>();
		targets.addAll(mapping.getListOfTargetMapping());
	
		ColumnReference cref = mapping.getColumnRef();
		if (cref == null) throw new NullPointerException("Column reference in an element mapping is NULL");
		
		ColumnDefinition col = table.getColumn(cref);
		if (col == null) throw new NullPointerException("A column declaration is undefined (name='" + cref.getColumnIdRef() + "')");
		
		mapping.getListOfTargetMapping().clear();
		if (DOSE.equals(col.getListOfColumnType().get(0))) {
			for (TargetMapping target : targets) {
				if (target == null) continue;
				int size = target.getListOfMap().size();
				if (size > 0) remapDoseTarget(table, mapping, target.getListOfMap());
				else throw new UnsupportedOperationException("Unrecognised dosing target cardinality (name='" + cref.getColumnIdRef() + "')");
			}
		} 
		else throw new UnsupportedOperationException("Unrecognised target mapping context on a data column (name='" + cref.getColumnIdRef() + "')");
	}
	
	private void remapDoseTarget(TabularDataset table, ColumnMapping mapping, List<MapType> maps) {
		
		if (mapping == null || maps == null || table == null) return;
		if (maps.isEmpty()) return;
		
		Piecewise pw = mapping.createPiecewise();
		
		for (MapType map : maps) {
			if (map == null) throw new NullPointerException("A dose target mapping record is NULL");
			Integer adm = map.getAdmNumber();
			String dataSymbol = map.getDataSymbol();
			
			if (adm == null) throw new NullPointerException("The administration number is not specified");
			
			CommonVariableDefinition v = getTranslatedElement(adm);
			if (v == null) throw new IllegalStateException("Unable to determine model element from the translated PK macro");
			
			if (dataSymbol != null) remapDoseTargetWithAdministrationAndCompartment(table, v, pw, dataSymbol);
			else remapDoseTargetWithAdministrationOnly(table, v, pw); 
		}
	}

	private void remapDoseTargetWithAdministrationAndCompartment(TabularDataset table, CommonVariableDefinition target, Piecewise pw, String dataSymbol) {
		if (table == null || target == null || pw == null || dataSymbol == null) return;
		ColumnDefinition adm_col = table.getColumn(ADM);
		ColumnDefinition dose_col = table.getColumn(DOSE);
		
		if (dose_col == null) throw new IllegalStateException("The dose column is not specified in the model external dataset declaration");
		if (adm_col == null) throw new IllegalStateException("The administration column (ADM) is not specified in the model external dataset declaration");
		
		Piece piece = addPiece(pw, true);
		piece.setValue(symbolRef(target, accessor));
		
		LogicBinOp ADM_clause = logicalBinaryOperator("eq", new Object[] {adm_col, Integer.valueOf(dataSymbol)}, accessor);
		LogicBinOp DOSE_clause = logicalBinaryOperator("gt", new Object[] {dose_col, 0}, accessor);
		LogicBinOp lbop = logicalBinaryOperator("and", new Object[] {ADM_clause, DOSE_clause}, accessor);
		piece.getCondition().setLogicBinop(lbop);
	}

	private void remapDoseTargetWithAdministrationOnly(TabularDataset table, CommonVariableDefinition target, Piecewise pw) {
		if (table == null || target == null || pw == null) return;
		ColumnDefinition dose_col = table.getColumn(DOSE);
		
		if (dose_col == null) throw new IllegalStateException("The dose column is not specified in the model external dataset declaration");
		
		Piece piece = addPiece(pw, true);
		piece.setValue(symbolRef(target, accessor));
		
		LogicBinOp DOSE_clause = logicalBinaryOperator("gt", new Object[] {dose_col, 0}, accessor);
		piece.getCondition().setLogicBinop(DOSE_clause);
	}

	/**
	 * Translate regressor fields to covariate.<br/>
	 * Convenience methods for non-Monolix language, i.e. everything else in the planet.
	 */
	private void resetRegressorToCovariate() {
		if (!resetRegToCov) return;
		
		TrialDesign td = dom.getTrialDesign();
		if (td == null) return;
		
		List<ExternalDataSet> exds = td.getListOfExternalDataSet();
		if (exds == null) return;
		if (exds.isEmpty()) return;
		
		for (ExternalDataSet exd : exds) {
			if (exd == null) continue;
			
			DataSet ds = exd.getDataSet();
			if (ds == null) continue;
			
			HeaderColumnsDefinition definition = ds.getDefinition();
			if (definition == null) continue;
			
			List<ColumnDefinition> cols = definition.getListOfColumn();
			if (cols == null) continue;
			for (ColumnDefinition col : cols) {
				if (col == null) continue;
				
				ColumnType type = col.getListOfColumnType().get(0);
				if (type == null) continue;
				
				if (type.equals(ColumnType.REG)) {
					col.setColumnType(ColumnType.COVARIATE);
					hasResetColumnUsageRegToCov = true;
				}
			}
		}
	}

	/**
	 * Set the accessor handle to a model.
	 * @param a Accessor handle
	 */
	public void setAccessor(Accessor a) { accessor = a; }
	
	@Override
	public void setAddPlottingBlock(boolean decision) { }

	@Override
	public void setCurrentStructuralBlock(StructuralBlock sb) {
		if (sb == null) throw new NullPointerException("Structural block is NULL.");
		currentSb = sb;
	}

	@Override
	public void setDeactivateIdFactory(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setDom(PharmML dom) { this.dom = dom; }

	@Override
	public void setDuplicateVariablesPermitted(boolean decision) { throw new UnsupportedOperationException(); }

	/**
	 * Instruct a converter to actively
	 * @param decision
	 */
	public void setEchoException(boolean decision) { is_echo_exception = decision; }

	/**
	 * Throw exception if input model has a continuous covariate.
	 * Otherwise the PFIM converter just returns an exception. 
	 * @param decision Decision
	 */
	public void setExceptionWithContinuousCovariate(boolean decision) { exceptionWithContinuousCovariate = decision; }
	
	@Override
	public void setFilterReservedWords(boolean decision) { filter_reserved_words = decision; }
	
	@Override
	public void setIndexFromZero(boolean decision) { throw new UnsupportedOperationException(); }
	
	@Override
	public void setIsolateConditionalDoseVariable(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setIsolateDoseTimingVariable(boolean decision) { isolate_dt = decision; }

	@Override
	public void setLexOnly(boolean decision) { lexOnly = decision; }

	@Override
	public void setLoadOnly(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setMixedEffect(boolean value) { throw new UnsupportedOperationException(); }
	
	/**
	 * Set the model filename if not explicitly set by filepath on the performConvert() method.
	 * @param model_filename_ Path/Name to the inputted model file.
	 */
	public void setModelFilename(String model_filename_) { model_filename = model_filename_; }

	/**
	 * Set the output directory for the converter instance.
	 * @param dir Output Directory
	 */
	public void setOutputDirectory(File dir) {
		boolean done = false;
		if (dir == null) throw new NullPointerException("Output directory is NULL.");
		if (dir.isDirectory()) {
			if (dir.canRead() && dir.canWrite()) {
				output_dir = dir.getAbsolutePath();
				done = true;
			}
		}

		if (!done) throw new IllegalStateException("Unable to assign output directory path.");
	}
	
	@Override
	public void setParser(IParser parser_) {
		if (parser_ == null) throw new NullPointerException("The parser is null.");
		else
			parser = (Parser) parser_;
	}
	
	@Override
	public void setPermitEmptyTrialDesignBlock(boolean decision) { permitEmptyTrialDesignBlock = decision; }
	
	@Override
	public void setRemoveIllegalCharacters(boolean decision) { remove_illegal_chars = decision; }

	/**
	 * Set the run identifier for output file stem.
	 * @param run_id_ Run Identifier
	 */
	public void setRunId(String run_id_) { run_id = run_id_; }
	
	@Override
	public void setSaveRenamedSymbolList(boolean decision) { save_renamed_symbol_list = true; }

	@Override
	public void setSaveSimulationOutput(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setScriptDefinition(ScriptDefinition sd_) { sd = sd_; }
	
	@Override
	public void setSortParameterModel(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setSortParameterModelByClustering(boolean decision) { throw new UnsupportedOperationException(); }
	
	@Override
	public void setSortParameterModelByContext(boolean decision) { throw new UnsupportedOperationException(); }
	
	@Override
	public void setSortStructuralModel(boolean decision) { sort_structural_model = decision;	}

	@Override
	public void setSortStructuralModelByClustering(boolean decision) {}

	@Override
	public void setSortVariabilityLevels(boolean decision) { throw new UnsupportedOperationException(); }

	@Override
	public void setTerminateWithInvalidXML(boolean decision) { throw new UnsupportedOperationException(); }

    @Override
	public void setTranslate(boolean decision) { translate_macros = decision; }
    
    @Override
	public void setTranslator(Translator tr_) { if (tr_ != null) tr = tr_; }

    @Override
	public void setTreeMaker(TreeMaker tm) { if (tm != null) this.tm = tm; }

    @Override
	public void setUseGlobalConditionalDoseVariable(boolean decision) { throw new UnsupportedOperationException(); }

    @Override
	public void setUsePiecewiseAsEvents(boolean decision) { usePiecewiseAsEvents = decision; }

    @Override
	public void setValidateXML(boolean decision) { throw new UnsupportedOperationException(); }

    private void sortElementOrdering() throws NullPointerException, IOException { 
		if (sort_structural_model) sortStructuralBlock(getStrucuturalBlock());
	}

    private void sortStructuralBlock(StructuralBlock sb) throws NullPointerException, IOException {
		if (sb == null) return;

		List<PharmMLRootType> list = sb.getListOfDeclarations();
		List<DependencyRef> refs = new ArrayList<DependencyRef>();
		Map<Object, DependencyRef> dep_map = new HashMap<Object, DependencyRef>();
		for (PharmMLRootType o : list) {
			DependencyRef ref = new DependencyRef(o);
			BinaryTree bt = tm.newInstance(o);
			
			List<NestedTreeRef> ntrefs = tm.getNestedTrees();
			List<NestedTreeRef> local_ntrefs = new ArrayList<NestedTreeRef>();
			local_ntrefs.addAll(ntrefs);
			
			addDependency(ref, bt);
			for (NestedTreeRef ntref : local_ntrefs) addDependency(ref, ntref.bt);
			
			refs.add(ref);
			dep_map.put(o, ref);
		}
		
		for (int i = 0; i < refs.size(); i++) {
			DependencyRef ref = refs.get(i);
			if (ref.hasDependendsUpon()) {
				for (int j = 0; j < ref.getDependsUpon().size(); j++) {
					PharmMLElement depends_upon = ref.getDependsUpon().get(j);			
					DependencyRef other_ref = dep_map.get(depends_upon);
					if (other_ref != null) {
						if (other_ref.hasDependendsUpon()) {
							for (int k = 0; k < other_ref.getDependsUpon().size(); k++) {
								ref.addDependency(other_ref.getDependsUpon().get(k));
							}
						}
					}
				}
			}
		}
		List<PharmMLElement> elements_under_consideration = createElementsUnderConsideration(refs);
		
		if (elements_under_consideration.isEmpty()) return;
		updateDependencyContext(elements_under_consideration, refs);

		DependencyGraph g = new DependencyGraph(refs);

		g.createVertices();
		int edgeCount = g.createEdges();
		if (edgeCount == 0) return;

		g.sort();

		List<PharmMLElement> ordered_variables = g.getSortedElements();
		if (ordered_variables.isEmpty()) return;

		sb.setOrderedVariableList(ordered_variables);
	}
    
	/**
	 * Function called to translate PK macros.<br/>
	 * Override to change application logic.
	 */
	private void translatePKMacros() {
		if (!translate_macros) return;
		
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		
		List<StructuralModel> sms_with_macros = new ArrayList<StructuralModel>(); 
		List<StructuralModel> sms = def.getListOfStructuralModel();
		
		Accessor a = new Accessor(dom);
		for (StructuralModel sm : sms) {
			if (sm == null) continue;
			
			PKMacroList macros = a.getPKMacros(sm);
			if (macros == null) continue;
			if (macros.getListOfMacro().isEmpty()) continue;
			sms_with_macros.add(sm);
		}
		
		if (sms_with_macros.size() == 0) return;
		
		if (tr == null) throw new NullPointerException("The Macro translator is NULL.");
		
		for (StructuralModel sm : sms_with_macros) {
			try {
				MacroOutput output = tr.translate(sm, target_level, accessor.getIndependentVariable());
				StructuralModel sm_translated = output.getStructuralModel();
				if(!replace(dom, sm, sm_translated))
					throw new IllegalStateException("PK macros translation failed (blkId='" + sm.getBlkId() + "')");
				
				// Register the macro input and outputs so that they can be associated 
				// with the translated structural block later on.
				PKMacroList macros = accessor.getPKMacros(sm);
				
				if (!macro_output_map.containsKey(sm_translated)) macro_output_map.put(sm_translated, output);
				if (!macro_input_map.containsKey(sm_translated)) macro_input_map.put(sm_translated, macros.getListOfMacro()); 
			} catch (InvalidMacroException e) {
				throw new UnsupportedOperationException(e);
			}
		}
		
		updatePKMacroDataMappings();
	}

	@Override
	public void updateNestedTrees() {
		for (NestedTreeRef ref : tm.getNestedTrees()) addStatement(ref);
		tm.getNestedTrees().clear();
	}

	private void updatePKMacroDataMappings() {
		TrialDesign td = dom.getTrialDesign();
		if (td == null) return;
		
		List<ExternalDataSet> exds = td.getListOfExternalDataSet();
		if (exds == null) return;
		
		ModelDefinition def = dom.getModelDefinition();
		if (def == null) return;
		getAccessor();
		
		for (ExternalDataSet exd : exds) {
			if (exd == null) continue;

			TabularDataset table = new TabularDataset(exd.getDataSet(), this);
			for (PharmMLRootType element : exd.getListOfColumnMappingOrColumnTransformationOrMultipleDVMapping()) {
				if (!isColumnMapping(element)) continue;
				remapColumnPostPKMacroTranslation((ColumnMapping) element, table);
			}
		}
	}

	@Override
	public boolean useCachedDependencyList() { throw new UnsupportedOperationException(); }

	@Override
	public List<Symbol> getContinuousOutputs() {
		List<Symbol> outputs = new ArrayList<Symbol>();
		List<ObservationBlock> obs = getObservationBlocks();
		
		for (ObservationBlock ob : obs) {
			if (ob == null) continue;
			ObservationModel model = ob.getModel();
			if (model == null) continue;
			ContinuousObservationModel com = model.getContinuousData();
			if (com == null) continue; 
			ObservationError oe = com.getObservationError();
			if (oe == null) continue;
			if (isStructuredError(oe)) {
				StructuredObsError soe = (StructuredObsError) oe;
				PharmMLRootType element = accessor.fetchElement(soe.getOutput());
				if (isSymbol(element)) outputs.add((Symbol) element); 
			}
		}
		
		return outputs;
	}
}

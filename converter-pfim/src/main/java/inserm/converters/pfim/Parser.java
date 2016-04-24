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

import static crx.converter.engine.BaseModelWriter.vector;
import static crx.converter.engine.PharmMLTypeChecker.isCommonParameter;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isVector;
import static inserm.converters.pfim.GenericOption.ALPHA;
import static inserm.converters.pfim.GenericOption.ATOL;
import static inserm.converters.pfim.GenericOption.BETA_COVARIATE;
import static inserm.converters.pfim.GenericOption.BETA_COVARIATE_OCCASSION;
import static inserm.converters.pfim.GenericOption.COMPUTE_NNI;
import static inserm.converters.pfim.GenericOption.COMPUTE_NNI_EQUIVALENCE;
import static inserm.converters.pfim.GenericOption.COMPUTE_POWER;
import static inserm.converters.pfim.GenericOption.COMPUTE_POWER_EQUIVALENCE;
import static inserm.converters.pfim.GenericOption.COVARIATE_OCCASSION_CATEGORIES;
import static inserm.converters.pfim.GenericOption.COVARIATE_OCCASSION_NAMES;
import static inserm.converters.pfim.GenericOption.COVARIATE_OCCASSION_PROPORTIONS;
import static inserm.converters.pfim.GenericOption.COVARIATE_OCCASSION_SEQUENCE;
import static inserm.converters.pfim.GenericOption.COVARIATE_PROPORTIONS;
import static inserm.converters.pfim.GenericOption.DOSE_IDENTICAL;
import static inserm.converters.pfim.GenericOption.DOSE_VECTOR;
import static inserm.converters.pfim.GenericOption.EQUIVALENCE_INTERVAL_START;
import static inserm.converters.pfim.GenericOption.EQUIVALENCE_INTERVAL_STOP;
import static inserm.converters.pfim.GenericOption.GAMMA;
import static inserm.converters.pfim.GenericOption.GIVEN_POWER;
import static inserm.converters.pfim.GenericOption.IDENTICAL_COND_4_EACH_DESIGN;
import static inserm.converters.pfim.GenericOption.IDENTICAL_TIMES;
import static inserm.converters.pfim.GenericOption.NR;
import static inserm.converters.pfim.GenericOption.NSUBJECTS;
import static inserm.converters.pfim.GenericOption.N_OCC;
import static inserm.converters.pfim.GenericOption.OMEGA;
import static inserm.converters.pfim.GenericOption.PARAMETER_ASSOCIATED;
import static inserm.converters.pfim.GenericOption.PARAMS_ASSOCIATED_WITH_OCCASSION;
import static inserm.converters.pfim.GenericOption.PROPORTIONS;
import static inserm.converters.pfim.GenericOption.PROTA;
import static inserm.converters.pfim.GenericOption.RTOL;
import static inserm.converters.pfim.GenericOption.SIG_INTERA;
import static inserm.converters.pfim.GenericOption.SIG_SLOPEA;
import static inserm.converters.pfim.GenericOption.TRAND;
import static inserm.converters.pfim.GenericOption.USING_COVARIATE;
import static inserm.converters.pfim.GenericOption.USING_COVARIATE_OCCASSION;
import static inserm.converters.pfim.GraphOption.GRAPH_LOGICAL;
import static inserm.converters.pfim.GraphOption.GRAPH_SENSITIVITY;
import static inserm.converters.pfim.GraphOption.SAMPLING_TIME_LOWER;
import static inserm.converters.pfim.GraphOption.SAMPLING_TIME_UPPER;
import static inserm.converters.pfim.GraphOption.STANDARD_GRAPHIC;
import static inserm.converters.pfim.GraphOption.X_AXES_NAMES;
import static inserm.converters.pfim.GraphOption.Y_AXES_NAMES;
import static inserm.converters.pfim.GraphOption.Y_AXES_RANGE;
import static inserm.converters.pfim.OptimisationAlgorithm.SIMPLEX;
import static inserm.converters.pfim.SimplexOption.DELTA_TIME;
import static inserm.converters.pfim.SimplexOption.MAX_ITERATIONS;
import static inserm.converters.pfim.SimplexOption.OPTIMISATION_OF_PROPORTIONS_OF_SUBJECTS;
import static inserm.converters.pfim.SimplexOption.PRINT_ITERATIONS;
import static inserm.converters.pfim.SimplexOption.RELATIVE_CONVERGENCE_TOLERANCE;
import static inserm.converters.pfim.SimplexOption.SAMPLING_TIME_LOWER_A;
import static inserm.converters.pfim.SimplexOption.SAMPLING_TIME_LOWER_B;
import static inserm.converters.pfim.SimplexOption.SAMPLING_TIME_UPPER_A;
import static inserm.converters.pfim.SimplexOption.SAMPLING_TIME_UPPER_B;
import static inserm.converters.pfim.SimplexOption.SIMPLEX_PARAMETER;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.engine.Accessor;
import crx.converter.engine.common.CategoryProportions;
import crx.converter.engine.common.CovariateParameterRef;
import crx.converter.engine.parts.CategoricalCovariateRef;
import crx.converter.engine.parts.CovariateBlockImpl;
import crx.converter.spi.blocks.CovariateBlock;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.steps.EstimationStep;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.BooleanValue;
import eu.ddmore.libpharmml.dom.commontypes.FalseBoolean;
import eu.ddmore.libpharmml.dom.commontypes.IntValue;
import eu.ddmore.libpharmml.dom.commontypes.Matrix;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLElement;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.Symbol;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.TrueBoolean;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.commontypes.Vector;
import eu.ddmore.libpharmml.dom.commontypes.VectorElements;
import eu.ddmore.libpharmml.dom.commontypes.VectorValue;
import eu.ddmore.libpharmml.dom.modeldefn.CategoricalCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.Category;
import eu.ddmore.libpharmml.dom.modeldefn.CommonParameter;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modellingsteps.Algorithm;
import eu.ddmore.libpharmml.dom.modellingsteps.EstimationOpType;
import eu.ddmore.libpharmml.dom.modellingsteps.EstimationOperation;
import eu.ddmore.libpharmml.dom.modellingsteps.OperationProperty;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;

public class Parser extends crx.converters.r.Parser {
	private static class BetaCovariate {
		public CovariateDefinition cov = null;
		public List<VectorValue> values = null;
		
		public BetaCovariate(CovariateDefinition cov_, Rhs rhs) {
			if (cov_ == null) throw new NullPointerException("Covariate is NULL");
			if (rhs == null) throw new NullPointerException("Rhs is NULL");
			
			cov = cov_;
			
			Object content = rhs.getContent();
			if (!isVector(content)) throw new IllegalStateException("Beta covariate parameters not stored in a vector.");
			
			Vector v = (Vector) content;
			values = v.getVectorElements().getListOfElements();
		}
	}
	
	private static class OccasionSequences {
		public List<Rhs> sequences = new ArrayList<Rhs>();
		
		public boolean addSequence(Rhs rhs) {
			if (rhs == null) return false;
			
			if (!sequences.contains(rhs)) {
				sequences.add(rhs);
				return true;
			}
			
			return false;
		}
	}
	
	private static class ParametersAssociatedWithOccassion {
		public List<CommonParameter> params = new ArrayList<CommonParameter>();
		
		public boolean addParameter(CommonParameter p) {
			if (p == null) return false;
			
			if (!params.contains(p)) {
				params.add(p);
				return true;
			}
			
			return false;
		}
	}
	
	private static final String CURRENT_WORKING_DIR = "CURRENT_WORKING_DIR";
	private static final String PFIM_PROGRAM_DIR = "PFIM_PROGRAM_DIR";
	private static final String pfimProjectFilename = "PFIM";
	private static String pfimStdinFilename = "stdin";
	private static String pfimStdoutFilename = "Stdout";
	private static final String PREFERRED_SEPERATOR = "/";
	
	private OptimisationAlgorithm algo = OptimisationAlgorithm.UNSPECIFIED;
	private Double alpha = 0.01;
	private Double atol = 1E-08;
	private List<BetaCovariate> beta_covariates = new ArrayList<BetaCovariate>();
	private List<CategoryProportions> cat_props = new ArrayList<CategoryProportions>();
	private boolean computeNNI = false;
	private boolean computeNNIEquivalence = false;
	private Boolean computePower = false;
	private Boolean computePowerEquivalence = false;
	private Vector covOccassionNames = null;
	private FIMOption currentFIMOption = FIMOption.Population;
	private Double deltaTime = 0.0;
	private boolean diagonalFIM = false;
	private Rhs doseVector = null;
	private Rhs equivalence_interval_start = null;
	private Rhs equivalence_interval_stop = null;
	private Matrix gamma = null;
	private Double givenPower = 0.0;
	private Boolean graphLogical = false;
	private Boolean graphLogicalSensitivity = false;
	private String hmax = "Inf";
	private boolean identicalInitialConditionsInEachDesign = false;
	private Boolean identicalTimes = false;
	private PharmMLElement initialConditionTime = new IntValue(0);
	private Vector lowerA = null, lowerB = null;
	private Vector lowerSamplingTimes = null;
	private Integer maximumIterations = 5000;
	private Integer numberOfOccassions = 0;
	private Integer numberOfResponses = 0;
	private Integer numberOfSubjects = 0;
	private Map<String, Vector> occ_beta_covariate_map = new HashMap<String, Vector>();
	private Map<String, ParametersAssociatedWithOccassion> occ_parameter_map = new HashMap<String, ParametersAssociatedWithOccassion>();
	private Map<String, OccasionSequences> occ_seq_map = new HashMap<String, OccasionSequences>();
	private List<String> occassionCategories = new ArrayList<String>();
	private Map<String, Vector> occassionProportions_map = new HashMap<String, Vector>();
	private Matrix omega = null;
	private Boolean optimisationOfProportionsOfSubjects = false;
	private String outputFIMFilename = "";
	private List<String> pfimProjectTemplate = new ArrayList<String>();
	private String previousFIM = "";
	private Boolean printIterations = true;
	private String programDirectory = ".";
	private List<Object> protA = new ArrayList<Object>();
	private Double relativeConvergenceTolerance = 1E-6;
	private Double rtol = 1E-08;
	private Double sigInterA = 0.0;
	private Double sigSlopeA = 0.0;
	private Double simplexParameter = 0.0;
	private Boolean standardGraphic = true;
	private TreeMaker tm = null;
	private Integer Trand = -1;
	private Vector upperA = null, upperB = null;
	private Vector upperSamplingTimes = null;
	private boolean usingCovariateModel = false; 
	private boolean usingCovariateOccassionModel = false;
	private Boolean usingIdenticalDose = false;
	private boolean usingSubjects = false;
	private boolean writingAnalyticalModel = false;
	private boolean writingModelFunctionScopedVariable = false;
	private boolean writtenSTDIN = false;
	private boolean wrotePFIM_R = false;
	private Vector xAxesNames = null;
	private Vector yAxesNames = null;
	private Vector yAxesRange = null;
	
	public Parser() throws IOException {
		super();
		loadPFIMTemplate();
		setUseCompiledStructuralModel(false);
		
		// PFIM Specific symbols and settings.
		param_model_symbol = "p";
		state_vector_symbol = "X";
		programDirectory = props.getProperty("pfimProgramDirectory");
		
		// Read directory path via the shell, ignoring the relative path
		// specified in the config file.
		if (System.getProperty("PFIM_PROG_DIR") != null) programDirectory = System.getProperty("PFIM_PROG_DIR");
		setDiagonalMatrixViaList(true);
	}
	
	@Override
	protected String doCategoricalCovariateRef(CategoricalCovariateRef cref) {
		String name = cref.cov.getSymbId();
		CategoricalCovariate cc = cref.cov.getCategorical();
		
		String format = "\"%s\"";
		StringBuffer categories = new StringBuffer();
		
		int i = 0;
		for (Category cat : cc.getListOfCategory()) {
			if (i > 0) categories.append(",");
			categories.append(parse(cat, lexer.getStatement(cat)).trim());
			i++;
		}
		
		format = "%s=c(%s)";
		return String.format(format, name, categories); 
	}
	
	@Override
	protected String doCovariateParameterRef(CovariateParameterRef ref) {
		String cname = z.get(ref.getCovariate());
		
		List<String> vars = new ArrayList<String>(); 
		for(CommonParameter p : ref.getParameters()) {
			if (p == null) continue;
			vars.add(getSymbol(new StringValue(z.get(p))));
		}
		
		 
		int i = 0;
		StringBuffer stmt = new StringBuffer(cname + "=c(");
		for (String var :  vars) {
			if (i > 0) stmt.append(",");
			stmt.append(var);
			i++;
		}
		stmt.append(")");
		
		return stmt.toString();
	}
	
	private void doPFIMLocalAssignment(PrintWriter fout, Symbol v) {
		if (fout  == null || v == null) return;
		
		String symbol = getSymbol(v);
		if (isLocalVariable(v)) symbol = state_vector_symbol;
		String format = "%s <- %s\n";
		
		fout.write(String.format(format, z.get(v.getSymbId()), symbol));
	}
	
	@Override
	protected String doSymbolRef(SymbolRef s) {
		// X symbol reserved symbol as used by PFIM.
		if (writingAnalyticalModel) {
			Accessor a = lexer.getAccessor();
			PharmMLRootType element = a.fetchElement(s);
			if (isLocalVariable(element)) return state_vector_symbol;
		}
		
		if (writingModelFunctionScopedVariable) return z.get(s.getSymbIdRef());
		else return super.doSymbolRef(s);
	}
	
	private String getBetaCovariateDecl(BetaCovariate bcov) {
		TreeMaker tm = lexer.getTreeMaker();
		boolean logTransform = false;
		if (Trand == 2) logTransform = true;
		
		String format = "log(%s)";
		StringBuffer stmt = new StringBuffer("c(");
		int i = 0;
		for (VectorValue v : bcov.values) {
			if (i > 0) stmt.append(",");
			String value = stripOuterBrackets(parse(bcov, tm.newInstance(v)).trim()); 
			if (logTransform) value = String.format(format, value);
			stmt.append(value);
			i++;
		}
		stmt.append(")");
		
		format = "%s=list(%s)";
		return String.format(format, z.get(bcov.cov), stmt);
	}
	
	@Override
	protected String getModelFunctionName(StructuralBlock sb) { return "model"; }
	
	private String getOccassionCategories(String name, Vector values) {
		if (name == null || values == null) return null;
		
		TreeMaker tm = lexer.getTreeMaker();
		BinaryTree bt = tm.newInstance(values);
		lexer.updateNestedTrees();
		
		String stmt = stripOuterBrackets(parse(values, bt));
		String format = "%s=%s";
		return String.format(format, name, stmt);
	}
	
	// Try to use the canonical path so that the 'R' will find a target script no matter where things
	// are run.
	private String getPFIMProjectFilepath() throws IOException {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + pfimProjectFilename + "." + script_file_suffix;
	}
	
	private String getStdinFilepath() {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + pfimStdinFilename + "." + script_file_suffix;
	}
	
	private void loadPFIMTemplate() {
		try {
			InputStream resource = getClass().getResourceAsStream("PFIM.r");
			if (resource == null) throw new NullPointerException("Unable to load PFIM template file.");
			
			BufferedReader in= new BufferedReader(new InputStreamReader(resource));
			String line = null; 
			while ((line = in.readLine()) != null) pfimProjectTemplate.add(line);
			
			in.close();
			in = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void readAlgorithm(EstimationOperation op) {
		if (op == null) return;
		
		Algorithm algorithm = op.getAlgorithm();
		if (algorithm == null) return;
		
		String definition = algorithm.getDefinition();
		if (definition == null) return;
		
		if(!OptimisationAlgorithm.contains(definition)) 
			throw new UnsupportedOperationException("The algorithm (definition='" +  definition + "') is not supported by PFIM.");
		
		algo = OptimisationAlgorithm.fromValue(definition);
		if (algo.equals(SIMPLEX)) readSimplexOptions(op);
	}
	
	private void readEstimationOperations() {
		EstimationStep est = lexer.getEstimationStep();
		EstimationOperation [] ops = est.getOperations();
		
		// Explicitly look for a FIM operation.
		List<EstimationOperation> operation_list = new ArrayList<EstimationOperation>();
		for (EstimationOperation op : ops) {
			if (op == null) continue;
			String type = op.getOpType();
			
			if (type == null) throw new NullPointerException("Estimation operation is NULL.");
			
			boolean opFound = false;
			for (EstimationOpType op_typedef : EstimationOpType.values()) {
				if (op_typedef.value().equals(type)) {
					operation_list.add(op);
					opFound = true;
				}
			}
			
			if (!opFound) throw new UnsupportedOperationException("Estimation operation not recognised (type='" + type + "')");
			
			readAlgorithm(op);
		}
		
		if (operation_list.isEmpty()) throw new IllegalStateException("No operations specified for the PFIM converter.");
		
		for (EstimationOperation operation : operation_list) {
			String type = operation.getOpType();
			if (type.equalsIgnoreCase(EstimationOpType.EST_POP.value())) readPopulationEstimationSettings(operation);
			else if (type.equalsIgnoreCase(EstimationOpType.EST_FIM.value())) readFIMEstimationSettings(operation); 
		}
	}
	
	private void readFIMEstimationSettings(EstimationOperation fim_est) {
		List<OperationProperty> props = fim_est.getProperty();
		
		for (OperationProperty prop : props) {
			if (prop == null) continue;
			String name = prop.getName();
			if (name == null) continue;
			
			if (name.equals("type")) setFIMType(prop);
			else if (name.equalsIgnoreCase("calculateFullFIM")) setFullFIMOption(prop);
		}
	}
	
	private void readGeneralOperationProperties(List<OperationProperty> props) {
		if (props == null) return;
		if (props.isEmpty()) return;
		
		for (OperationProperty prop : props) {
			if (prop == null) continue;
			String name = prop.getName();
			
			if (!GenericOption.contains(name) && !GraphOption.contains(name)) continue;
			
			if (name.equals(ATOL.toString())) setAtol(prop);
			else if (name.equals(RTOL.toString())) setRtol(prop);
			else if (name.equals(N_OCC.toString())) setNumberOfOccassions(prop);
			else if (name.equals(NR.toString())) setNumberOfResponses(prop);
			else if (name.equals(TRAND.toString())) setRandomEffectModel(prop);
			else if (name.equals(GAMMA.toString())) setGamma(prop);
			else if (name.equals(OMEGA.toString())) setOmega(prop);
			else if (name.equals(SIG_INTERA.toString())) setSigInterA(prop);
			else if (name.equals(SIG_SLOPEA.toString())) setSigSlopeA(prop);
			else if (name.equals(PROTA.toString())) setProtA(prop);
			else if (name.equals(NSUBJECTS.toString())) setNumberOfSubjects(prop);
			else if (name.equals(PROPORTIONS.toString())) setUsingProportions(prop);
			else if (name.equals(USING_COVARIATE.toString())) setUsingCovariateModel(prop);
			else if (name.equals(DOSE_IDENTICAL.toString())) setDoseIdentical(prop);
			else if (name.equals(DOSE_VECTOR.toString())) setDoseVector(prop);
			else if (name.equals(IDENTICAL_COND_4_EACH_DESIGN.toString())) setIdenticalInitialConditionsInEachDesign(prop);
			else if (name.startsWith(COVARIATE_PROPORTIONS.toString())) setCategoryProportions(prop);
			else if (name.startsWith(PARAMETER_ASSOCIATED.toString())) setParameterAssociatedWithCovariate(prop);
			else if (name.startsWith(BETA_COVARIATE.toString())) setBetaCovariate(prop);
			else if (name.startsWith(COVARIATE_PROPORTIONS.toString())) setCategoryProportions(prop);
			else if (name.equals(USING_COVARIATE_OCCASSION.toString())) setUsingCovariateOccassionModel(prop);
			else if (name.equals(COVARIATE_OCCASSION_NAMES.toString())) setCovariateOccassionNames(prop);
			else if (name.startsWith(COVARIATE_OCCASSION_CATEGORIES.toString())) setCovariateOccassionCategories(prop);
			else if (name.startsWith(COVARIATE_OCCASSION_SEQUENCE.toString())) setCovariateOccassionSequence(prop);
			else if (name.startsWith(COVARIATE_OCCASSION_PROPORTIONS.toString())) setOccassionProportions(prop);	
			else if (name.startsWith(PARAMS_ASSOCIATED_WITH_OCCASSION.toString())) setParameterAssociatedWithOccassion(prop);
			else if (name.startsWith(BETA_COVARIATE_OCCASSION.toString())) setBetaCovariateOccassion(prop);
			else if (name.equals(ALPHA.toString())) setAlpha(prop);
			else if (name.equals(COMPUTE_POWER.toString())) setComputePower(prop);
			else if (name.equals(COMPUTE_NNI.toString())) setNNI(prop);
			else if (name.equals(EQUIVALENCE_INTERVAL_START.toString())) equivalence_interval_start = prop.getAssign();
			else if (name.equals(EQUIVALENCE_INTERVAL_STOP.toString())) equivalence_interval_stop = prop.getAssign();
			else if (name.equals(COMPUTE_POWER_EQUIVALENCE.toString())) setComputePowerEquivalence(prop);
			else if (name.equals(COMPUTE_NNI_EQUIVALENCE.toString())) setComputeNNIEquivalence(prop);
			else if (name.equals(GIVEN_POWER.toString())) setGivenPower(prop);
			else if (name.equals(IDENTICAL_TIMES.toString())) setIdenticalTimes(prop);
			else if (name.equals(GRAPH_LOGICAL.toString())) setGraphLogical(prop);
			else if (name.equals(GRAPH_SENSITIVITY.toString())) setGraphLogicalSensitivity(prop);
			else if (name.equals(X_AXES_NAMES.toString())) setXAxesNames(prop);
			else if (name.equals(Y_AXES_NAMES.toString())) setYAxesNames(prop);
			else if (name.equals(STANDARD_GRAPHIC.toString())) setStandardGraphic(prop);
			else if (name.equals(SAMPLING_TIME_LOWER.toString())) setUpperSamplingTime(prop);
			else if (name.equals(SAMPLING_TIME_UPPER.toString())) setLowerSamplingTime(prop);
			else if (name.equals(Y_AXES_RANGE.toString())) setYAxesRange(prop);
		}
	}
	
	private void readOperations() { if (lexer.hasEstimation()) readEstimationOperations(); }
	
	private void readPopulationEstimationSettings(EstimationOperation pop_est) {
		if (pop_est == null) return;
		readGeneralOperationProperties(pop_est.getProperty());
	}
	
	private void readSimplexOptions(EstimationOperation op) {
		for (OperationProperty prop : op.getProperty()) {
			if (prop == null) continue;
			String name = prop.getName();
			if (!SimplexOption.contains(name)) continue;
		
			if (name.equals(DELTA_TIME.toString())) setDeltaTime(prop);
			else if (name.equals(MAX_ITERATIONS.toString())) setMaximumIterations(prop);
			else if (name.equals(OPTIMISATION_OF_PROPORTIONS_OF_SUBJECTS.toString())) setOptimisationOfProportionsOfSubjects(prop);
			else if (name.equals(PRINT_ITERATIONS.toString())) setPrintIterations(prop);
			else if (name.equals(RELATIVE_CONVERGENCE_TOLERANCE.toString())) setRelativeConvergenceTolerance(prop);
			else if (name.equals(SAMPLING_TIME_LOWER_A.toString())) setLowerA(prop);
			else if (name.equals(SAMPLING_TIME_LOWER_B.toString())) setLowerB(prop);
			else if (name.equals(SAMPLING_TIME_UPPER_A.toString())) setUpperA(prop);
			else if (name.equals(SAMPLING_TIME_UPPER_B.toString())) setUpperB(prop);
			else if (name.equals(SIMPLEX_PARAMETER.toString())) setSimplexParameter(prop);		
		}
	}
	
	private void setAlpha(OperationProperty prop) { alpha = readDouble(prop); }
	
	private void setAtol(OperationProperty prop) { atol = readDouble(prop); }
	
	public void setBetaCovariate(OperationProperty prop) {
		if (prop == null) return;
		Accessor a = lexer.getAccessor();
		String name = prop.getName().split("_")[1];
		CovariateDefinition cov = (CovariateDefinition) a.fetchElement(name);
		beta_covariates.add(new BetaCovariate(cov, prop.getAssign()));
	}
	
	private void setBetaCovariateOccassion(OperationProperty prop) {
		String name = prop.getName().split("_")[1];
		Object content = prop.getAssign().getContent();
		if (!isVector(content)) return;
		
		occ_beta_covariate_map.put(name, (Vector) content); 
	}
	
	private void setCategoryProportions(OperationProperty prop) {
		if (prop == null) return;
		
		Accessor a = lexer.getAccessor();
		String name = prop.getName().split("_")[1];
		CovariateDefinition cov = (CovariateDefinition) a.fetchElement(name);
		cat_props.add(new CategoryProportions(cov, prop.getAssign()));
	}
	
	private void setComputeNNIEquivalence(OperationProperty prop) { computeNNIEquivalence = readBoolean(prop); }
	
	private void setComputePower(OperationProperty prop) { computePower = readBoolean(prop); }
	
	private void setComputePowerEquivalence(OperationProperty prop) { computePowerEquivalence = readBoolean(prop); }
	
	private void setCovariateOccassionCategories(OperationProperty prop) {
		if (prop == null) return;
		
		Object content = prop.getAssign().getContent();
		if (!isVector(content)) return;
		Vector values = (Vector) content;
		
		String name = prop.getName().split("_")[1];
		occassionCategories.add(getOccassionCategories(name, values));
	}
	
	public void setCovariateOccassionNames(OperationProperty prop) {
		if (prop == null) return;
		
		Object content = prop.getAssign().getContent();
		if (!isVector(content)) return;
		covOccassionNames = (Vector) content;
	}
	
	private void setCovariateOccassionSequence(OperationProperty prop) {
		if (prop == null) return;
		
		String name = prop.getName().split("_")[1];
		
		OccasionSequences oseq = null;
		if (occ_seq_map.containsKey(name)) oseq = occ_seq_map.get(name);
		else {
			oseq = new OccasionSequences();
			occ_seq_map.put(name, oseq);
		}
		
		oseq.addSequence(prop.getAssign());
	}
	
	private void setDeltaTime(OperationProperty prop) { deltaTime = readDouble(prop); }
	
	private void setDoseIdentical(OperationProperty prop) {
		usingIdenticalDose = readBoolean(prop);
	}
	
	private void setDoseVector(OperationProperty prop) {
		doseVector = prop.getAssign();	
	}
	
	private void setFIMType(OperationProperty prop) { currentFIMOption = FIMOption.valueOf(readString(prop, true)); }
	
	private void setFullFIMOption(OperationProperty prop) {
		diagonalFIM = readBoolean(prop); 
	}
	
	private void setGamma(OperationProperty prop) {
		if (prop == null) return;
		gamma = prop.getAssign().getMatrix();
	}
	
	private void setGivenPower(OperationProperty prop) { givenPower = readDouble(prop); }
	
	private void setGraphLogical(OperationProperty prop) { graphLogical = readBoolean(prop); }
	
	private void setGraphLogicalSensitivity(OperationProperty prop) { graphLogicalSensitivity = readBoolean(prop); }
	
	private void setIdenticalInitialConditionsInEachDesign(OperationProperty prop) { identicalInitialConditionsInEachDesign = readBoolean(prop); }
	
	private void setIdenticalTimes(OperationProperty prop) { identicalTimes = readBoolean(prop); }
	
	private void setLowerA(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) lowerA = (Vector) content;
	}
	
	private void setLowerB(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) lowerB = (Vector) content;
	}
	
	private void setLowerSamplingTime(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) lowerSamplingTimes = (Vector) content;
	}
	
	private void setMaximumIterations(OperationProperty prop) { maximumIterations = readInt(prop); }
	
	private void setNNI(OperationProperty prop) { computeNNI = readBoolean(prop); }
	
	private void setNumberOfOccassions(OperationProperty prop) { numberOfOccassions = readInt(prop); }
	
	private void setNumberOfResponses(OperationProperty prop) { numberOfResponses = readInt(prop); }
	
	private void setNumberOfSubjects(OperationProperty prop) {
		usingSubjects = true;
		numberOfSubjects = readInt(prop);
	}
	
	private void setOccassionProportions(OperationProperty prop) {
		if (prop == null) return;
		
		String name = prop.getName().split("_")[1];
		if (!occassionProportions_map.containsKey(name)) {
			Object content = prop.getAssign().getContent();
			if (isVector(content)) occassionProportions_map.put(name, (Vector) content);
		}
	}
	
	private void setOmega(OperationProperty prop) {
		if (prop == null) return;
		omega = prop.getAssign().getMatrix();
	}
	
	private void setOptimisationOfProportionsOfSubjects(OperationProperty prop) { optimisationOfProportionsOfSubjects = readBoolean(prop); }
	
	public void setOutputFIMFilename(String filename) {
		if (filename != null) outputFIMFilename = filename;
	}
	
	private void setParameterAssociatedWithCovariate(OperationProperty prop) {
		if (prop == null) return;
		
		CovariateBlock cb = lexer.getCovariateBlock();
		if (cb == null) return;
		
		String cov_name = prop.getName().split("_")[1];
		Vector refs = (Vector) prop.getAssign().getContent();
		if (refs == null) return;
		
		VectorElements ve = refs.getVectorElements();
		if (ve == null) return;
		
		for (VectorValue value : ve.getListOfElements()) {
			if (isSymbolReference(value)) {
				SymbolRef ref = (SymbolRef) value;
				cb.addParameterToCovariate(cov_name, ref.getSymbIdRef());
			}
		}
	}
	
	private void setParameterAssociatedWithOccassion(OperationProperty prop) {
		if (prop == null) return;
		
		String name = prop.getName().split("_")[1];
		Object content = prop.getAssign().getContent();
		if (!isVector(content)) return;
		
		Accessor a = lexer.getAccessor();
		Vector v = (Vector) content;
		VectorElements elements = v.getVectorElements();
		ParametersAssociatedWithOccassion pawo = new ParametersAssociatedWithOccassion();
		for (VectorValue value : elements.getListOfElements()) {
			if (!isSymbolReference(value)) continue;
			SymbolRef ref = (SymbolRef) value;
			PharmMLRootType modelElement = a.fetchElement(ref);
			if (isCommonParameter(modelElement)) pawo.addParameter((CommonParameter) modelElement);
		}
		
		occ_parameter_map.put(name, pawo);
	}
	
	private void setPrintIterations(OperationProperty prop) { printIterations = readBoolean(prop); }
	
	private void setProtA(OperationProperty prop) {
		protA.add(prop.getAssign().getContent()); 
	}
	
	private void setRandomEffectModel(OperationProperty prop) { Trand = readInt(prop); }
	
	private void setRelativeConvergenceTolerance(OperationProperty prop) { relativeConvergenceTolerance = readDouble(prop); }
	
	private void setRtol(OperationProperty prop) { rtol = readDouble(prop); }
	
	private void setSigInterA(OperationProperty prop) { sigInterA = readDouble(prop); }
	
	private void setSigSlopeA(OperationProperty prop) { sigSlopeA = readDouble(prop); }
	
	private void setSimplexParameter(OperationProperty prop) { simplexParameter = readDouble(prop); }
	
	private void setStandardGraphic(OperationProperty prop) { standardGraphic = readBoolean(prop); }
	
	private void setUpperA(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) upperA = (Vector) content;
	}
	
	private void setUpperB(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) upperB = (Vector) content;
	}
	
	private void setUpperSamplingTime(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) upperSamplingTimes = (Vector) content;
	}
	
	private void setUsingCovariateModel(OperationProperty prop) {
		usingCovariateModel = readBoolean(prop);
	}
	
	public void setUsingCovariateOccassionModel(OperationProperty prop) {
		if (prop == null) return;
		usingCovariateOccassionModel = readBoolean(prop);
	}
	
	private void setUsingProportions(OperationProperty prop) {
		boolean result = readBoolean(prop);
		if (result) usingSubjects = false;
	}
	
	private void setXAxesNames(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) xAxesNames = (Vector) content;
	}
	
	private void setYAxesNames(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) yAxesNames = (Vector) content;
	}
	
	private void setYAxesRange(OperationProperty prop) { 
		Object content = prop.getAssign().getContent();
		if (isVector(content)) yAxesRange = (Vector) content;
	}
	
	private void writeAlgorithmOption(PrintWriter fout) {
		String format = "\nalgo.option<-\"%s\" %s Optimisation Algorithm\n";
		fout.write(String.format(format, algo.toString(), comment_char));
	}
	
	private void writeAlpha(PrintWriter fout) {
		String format = "\nalpha<-%s %s Type one error alpha\n";
		fout.write(String.format(format, alpha, comment_char));
	}
	
	private void writeAnalyticalFunctionPrototype(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		String format = "form<-function(%s,%s,%s) {\n";
		String iv = z.get(lexer.getAccessor().getIndependentVariable());
		
		fout.write(String.format(format, iv, param_model_symbol, state_vector_symbol));
	}
	
	private void writeAnalyticalModelFunction(PrintWriter fout, StructuralBlock sb) {
		writingAnalyticalModel = true;
		writeAnalyticalFunctionPrototype(fout, sb);
		writeAnalyticalModelFunctionContent(fout, sb);
		writingAnalyticalModel = false;
	}
	
	private void writeAnalyticalModelFunctionContent(PrintWriter fout, StructuralBlock sb) {
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb != null) for (PopulationParameter p : pb.getParameters()) doPFIMLocalAssignment(fout, p);
		
		// Assuming only a single variable
		List<VariableDefinition> locals = sb.getLocalVariables();
		if (locals.isEmpty()) throw new IllegalStateException("Model variable unspecified.");
		
		writingModelFunctionScopedVariable = true;
		VariableDefinition v = locals.get(0);
		parse(v, lexer.getStatement(v), fout);
		writingModelFunctionScopedVariable = false;
		
		String format = "return(%s)\n";
		fout.write(String.format(format, z.get(v)));
		
		fout.write("}\n");
	}
	
	private void writeAnalyticalModelOption(PrintWriter fout) {
		writeSectionHeader(fout, "ANALYTICAL MODEL OPTION");
		writeIdenticalDose(fout);
		writeDose(fout);
		writeFunctionTimeIntervals(fout);
		writeNumericalDerivative(fout);
	}
	
	private void writeBetaCovariate(PrintWriter fout) {
		if (beta_covariates.isEmpty()) return;
		
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		for (BetaCovariate bcov : beta_covariates) {
			if (i > 0) stmt.append(",");
			stmt.append(getBetaCovariateDecl(bcov));
			i++;
		}

		String format = "\nbeta.covariate<-list(%s)\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writeBetaCovariateOccassion(PrintWriter fout) {
		String format = "log(%s)";
		boolean logTransform = false;
		if (Trand == 2) logTransform = true;
		
		StringBuffer outer_stmt = new StringBuffer();
		int i = 0;
		
		for (String name : occ_beta_covariate_map.keySet()) {
			if (i > 0) outer_stmt.append(",");
			
			Vector v = occ_beta_covariate_map.get(name);
		
			int j = 0;
			StringBuffer inner_stmt = new StringBuffer(name + "=list(c(");
			List<VectorValue> values = v.getVectorElements().getListOfElements();
			for (VectorValue value : values) {
				if (j > 0) inner_stmt.append(",");
				String stmt = stripOuterBrackets(parse(v, tm.newInstance(value)).trim());
				if (logTransform) stmt = String.format(format, stmt);
				inner_stmt.append(stmt);
				j++;
			}
			inner_stmt.append("))");
			
			outer_stmt.append(inner_stmt);
			i++;
		}
		
		format = "\nbeta.covariate_occ<-list( %s )\n";
		fout.write(String.format(format, outer_stmt));
	}
	
	private void writeBetaFixed(PrintWriter fout) {
		EstimationStep est = lexer.getEstimationStep();
		if (est == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		List<PopulationParameter> ps = pb.getParameters();
		if (ps.isEmpty()) return;
		
		List<Object> values = new ArrayList<Object>();
		for (PopulationParameter p : ps) {
			if (p == null) continue;
			if (est.isFixedParameter(p)) values.add(new TrueBoolean());
			else values.add(new FalseBoolean());
		}
		
		Vector list = vector(values, lexer.getAccessor());
		
		BinaryTree bt = tm.newInstance(list);
		lexer.updateNestedTrees();
		
		String stmt = stripOuterBrackets(parse(list, bt).trim());
		String format = "\nbeta.fixed<-%s %s Fixed Effect Parameter names\n";
		fout.write(String.format(format, stmt, comment_char));
	}
	
	@Override
	public void writeCategoricalEstimationBlock(PrintWriter fout, File output_dir, EstimationStep est) throws IOException {
		writeStdin();
	}
	
	private void writeComputeNNIEquivalence(PrintWriter fout) {
		String value = parse(computeNNIEquivalence, tm.newInstance(computeNNIEquivalence));
		String format = "\ncompute.nni_eq<-%s %s Compute No. subjects needed for given equivalence test.\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeComputePower(PrintWriter fout) {
		String value = parse(computePower, tm.newInstance(computePower)).trim();
		
		String format = "\ncompute.power<-%s %s Compute expected comparison test power.\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeComputePowerEquivalence(PrintWriter fout) {
		String value = parse(computePowerEquivalence, tm.newInstance(computePowerEquivalence));
		String format = "\ncompute.power_eq<-%s %s Compute expected power for equivalence test\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeConditionsInitialises(PrintWriter fout) {
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb == null) return;
		
		if (sb.isPlainFunction()) writeConditionsInitialisesForAnalyticalForm(fout);
	}
	
	private void writeConditionsInitialisesForAnalyticalForm(PrintWriter fout) {
		List<Object> values = new ArrayList<Object>();
		values.add(30);
		values.add(0.0);
		
		PharmMLElement defaultInitialConditions = vector(values, lexer.getAccessor());
		
		BinaryTree bt = tm.newInstance(defaultInitialConditions);
		lexer.updateNestedTrees();
		
		String stmt = stripOuterBrackets(parse(defaultInitialConditions, bt).trim());
		String format = "\ncondinit<-%s\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writeCovariateCategoryList(PrintWriter fout) {
		CovariateBlockImpl cb = (CovariateBlockImpl) lexer.getCovariateBlock();
		if (cb == null) return;
		
		List<CategoricalCovariateRef> ccovs = cb.getCategoricalCovariates();
		if (ccovs.isEmpty()) return;
		
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		for (CategoricalCovariateRef ccov : ccovs) {
			if (ccov == null) continue;
			if (i > 0) stmt.append(",");
			stmt.append(parse(ccov, lexer.getStatement(ccov)));
			i++;
		}
		
		String format = "\ncovariate.category<-list(%s)\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writeCovariateModel(PrintWriter fout) {
		writeSectionHeader(fout, "Covariate Model");
		writeUsingCovariates(fout);
		writeCovariateNameVector(fout);
		writeCovariateCategoryList(fout);
		writeCovariateProportions(fout);
		writeParameterAssociated(fout);
		writeBetaCovariate(fout);
	}
	
	// Assuming just a single covariate model with the PharmML.
	private void writeCovariateNameVector(PrintWriter fout) {
		List<CovariateBlock> cbs = lexer.getCovariateBlocks();
		if (cbs == null) return;
		if (cbs.isEmpty()) return;
		
		CovariateBlock cb = cbs.get(0);
		if (cb == null) return;
		
		List<String> names = cb.getCategoricalCovariateNames();
		if (names.isEmpty()) return;
		
		List<Object> onames = new ArrayList<Object>();
		for (String name : names) onames.add(new StringValue(name));
		
		Vector cov_names = vector(onames, lexer.getAccessor());
		
		BinaryTree bt = tm.newInstance(cov_names);
		lexer.updateNestedTrees();
		
		String stmt = parse(cov_names, bt).trim();
		
		String format = "\ncovariate.name<-list(%s) %s Vector of covariates\n";
		fout.write(String.format(format, stripOuterBrackets(stmt), comment_char));
	}
	
	private void writeCovariateOccassionCategories(PrintWriter fout) {
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		for (String category : occassionCategories) {
			if (category == null) continue;
			if (i > 0) stmt.append(",");
			stmt.append(category);
			i++;
		}
		
		String format = "\ncovariate_occ.category<-list( %s ) %s 1st category the reference.\n";
		fout.write(String.format(format, stmt, comment_char));
	}
	
	private void writeCovariateOccassionModel(PrintWriter fout) {
		String value = parse(usingCovariateOccassionModel, tm.newInstance(usingCovariateOccassionModel)).trim();
		
		String format = "\ncovariate_occ.model<-%s %s Add covariate 'occassion' model\n";
		fout.write(String.format(format, value, comment_char));;
	}
	
	private void writeCovariateOccassionNames(PrintWriter fout) {
		if (covOccassionNames == null) return;
		
		BinaryTree bt = tm.newInstance(covOccassionNames);
		lexer.updateNestedTrees();
		
		String stmt = stripOuterBrackets(parse(covOccassionNames, bt));
		String format = "\ncovariate_occ.name<-list( %s ) %s Covariates depending on occassion\n";
		fout.write(String.format(format, stmt, comment_char));
	}
	
	private void writeCovariateOccassionProportions(PrintWriter fout) {
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		
		String format = "%s=%s";
		for (String name : occassionProportions_map.keySet()) {
			if (i > 0) stmt.append(",");
			BinaryTree bt = tm.newInstance(occassionProportions_map.get(name));
			lexer.updateNestedTrees();
			stmt.append(String.format(format, name, parse(name, bt).trim()));
			i++;
		}
		
		format = "\ncovariate_occ.proportions<-list( %s )\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writeCovariateOccassionSequences(PrintWriter fout) {
		StringBuffer outer_stmt = new StringBuffer();
		int i = 0;
		for (String occ : occ_seq_map.keySet()) {
			if (i > 0) outer_stmt.append(",");
			OccasionSequences oseq = occ_seq_map.get(occ);
			StringBuffer inner_stmt = new StringBuffer(occ + "=list(");
			int j = 0;
			for (Rhs rhs : oseq.sequences) {
				if (j > 0) inner_stmt.append(",");
				
				BinaryTree bt = tm.newInstance(rhs);
				lexer.updateNestedTrees();
				inner_stmt.append(parse(rhs, bt).trim());
				
				j++;
			}
			inner_stmt.append(")");
			
			outer_stmt.append(inner_stmt);
			i++;
		}
		
		String format = "\ncovariate_occ.sequence<-list( %s )\n";
		fout.write(String.format(format, outer_stmt));
	}

	private void writeCovariateProportions(PrintWriter fout) {
		if (cat_props.isEmpty()) return;
		
		List<String> decls = new ArrayList<String>();
		for (CategoryProportions cat_prop : cat_props) decls.add(parse(cat_prop, tm.newInstance(cat_prop)).trim());
		
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		for (String decl : decls) {
			if (i > 0) stmt.append(",");
			stmt.append(decl);
			i++;
		}
		
		String format = "\ncovariate.proportions<-list(%s)\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writeCovariatesChangingWithOccassion(PrintWriter fout) {
		writeSectionHeader(fout, "Covariates changing with occassion");
		writeCovariateOccassionModel(fout);
		writeCovariateOccassionNames(fout);
		writeCovariateOccassionCategories(fout);
		writeCovariateOccassionSequences(fout);
		writeCovariateOccassionProportions(fout);
		writeOccassionParameters(fout);
		writeBetaCovariateOccassion(fout);
	}
	
	private void writeDeltaTime(PrintWriter fout) {
		if (deltaTime == null) return;

		String format = "\ndelta.time<-%s %s Minimum time between 2 sampling times\n";
		fout.write(String.format(format, deltaTime, comment_char));
	}
	
	private void writeDifferentialEquationOption(PrintWriter fout) {
		writeSectionHeader(fout, "DIFFERENTIAL EQUATION OPTION");
		writeInitialConditionTime(fout);
		writeInitialConditionsIdentical(fout);
		writeConditionsInitialises(fout);
		writeSolverSettings(fout);
	}
	
	private void writeDose(PrintWriter fout) {
		if (doseVector == null) return;
		
		BinaryTree bt = tm.newInstance(doseVector);
		lexer.updateNestedTrees();
		
		String value = parse(doseVector, bt).trim();
		String format = "\ndose<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeEffectsAndErrorModel(PrintWriter fout) {
		writeSectionHeader(fout, "FIXED/RANDOM EFFECTS AND ERROR MODEL");
		writeFixedEffectParameterNames(fout);
		writeFixedEffectParameterValues(fout);
		writeBetaFixed(fout);
		writeNumberOfOccassions(fout);
		writeRandomEffectModel(fout);
		writeOmega(fout);
		writeGamma(fout);
		writeResidualError(fout);
		writeProtA(fout);
		writeSubjects(fout);
		writeProportions(fout);
	}
	
	private void writeEquivalenceInterval(PrintWriter fout) {
		if (equivalence_interval_start == null || equivalence_interval_stop == null) return;
		
		Object ctx = new Object();
		String start = parse(ctx, tm.newInstance(equivalence_interval_start)).trim();
		String stop = parse(ctx, tm.newInstance(equivalence_interval_stop)).trim();
		
		String format = "\ninterval_eq<-c(%s,%s) %s Equivalence Interval\n";
		fout.write(String.format(format, start, stop, comment_char));
	}
	
	//@Override
	public void writeEstimationBlock(PrintWriter fout, File output_dir, EstimationStep est) throws IOException {
		writeStdin();
		
		if (fout == null) return;
		String format = "\n%s()\n";
		fout.write(String.format(format, pfimProjectFilename));
	}
	
	private void writeFIM(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "\nFIM<-\"%s\"\n";
		fout.write(String.format(format, currentFIMOption));
	}
	
	private void writeFimTypeOption(PrintWriter fout) {
		String option = "2";
		String comment = "Complete FIM";
		
		if (diagonalFIM) {
			option = "1";
			comment = "Diagonal FIM";
		}
		
		String format = "\noption<-%s %s %s\n";
		fout.write(String.format(format, option, comment_char, comment));
	}
	
	private void writeFixedEffectParameterNames(PrintWriter fout) {
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		List<PopulationParameter> ps = pb.getParameters();
		if (ps.isEmpty()) return;
		
		List<Object> names = new ArrayList<Object>();
		for (PopulationParameter p : ps) {
			if (p == null) continue; 
			names.add(new StringValue(p.getSymbId()));
		}
		
		Vector list = vector(names, lexer.getAccessor());
		
		TreeMaker tm = lexer.getTreeMaker();
		BinaryTree bt = tm.newInstance(list);
		lexer.updateNestedTrees();
		
		String stmt = stripOuterBrackets(parse(list, bt).trim());
		String format = "\nparameters<-%s %s Fixed Effect Parameter names\n";
		fout.write(String.format(format, stmt, comment_char));
	}
	
	private void writeFixedEffectParameterValues(PrintWriter fout) {
		ParameterBlock pm = lexer.getParameterBlock();
		EstimationStep est = lexer.getEstimationStep();
		
		if (pm == null || est == null) return;
		
		List<PopulationParameter> params = pm.getParameters();
		if (params.isEmpty()) return;
		
		List<String> values = new ArrayList<String>();
		
		for (PopulationParameter p : params) {
			if (p == null) throw new NullPointerException("A population parameter is NULL");
			
			ParameterEstimate pe = est.getParameterEstimate(p);
			if (pe == null) continue;
			
			BinaryTree bt = tm.newInstance(pe.getInitialEstimate());
			lexer.updateNestedTrees();
			values.add(parse(pe, bt).trim());
		}
		
		StringBuffer stmt = new StringBuffer("c(");
		int i = 0;
		for (String value : values) {
			if (i > 0) stmt.append(",");
			stmt.append(value);
			i++;
		}
		stmt.append(")");
		
		String format = "\nbeta<-%s %s Fixed effect parameter values\n";
		fout.write(String.format(format, stmt, comment_char));
	}
	
	private void writeFunctionTimeIntervals(PrintWriter fout) {
		String format = "\n%s Time intervals of each expression\n";
		fout.write(String.format(format, comment_char));
		
		format = "bound%s<-list(c(0,Inf))\n";
		
		char c = 'A'; // One statement per response.
		for (int i = 1; i <= numberOfResponses; i++, c++) fout.write(String.format(format, c));
	}
	
	private void writeGamma(PrintWriter fout) {
		if (gamma == null) return;
		
		BinaryTree bt = tm.newInstance(gamma);
		String value = parse(gamma, bt).trim();
		
		String format = "\ngamma<-%s %s Inter-occassion random effects\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeGivenPower(PrintWriter fout) {
		String format = "\ngiven.power<-%s\n";
		fout.write(String.format(format, givenPower));
	}
	
	private void writeGraphLogical(PrintWriter fout) {
		String value = stripOuterBrackets(parse(graphLogical, tm.newInstance(graphLogical)));
		String format = "\ngraph.logical<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphLogicalSensitivity(PrintWriter fout) {
		String value = stripOuterBrackets(parse(graphLogicalSensitivity, tm.newInstance(graphLogicalSensitivity)));
		String format = "\ngraphsensi.logical<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphLowerSamplingTimes(PrintWriter fout) {
		if (lowerSamplingTimes == null) return;

		BinaryTree bt = tm.newInstance(lowerSamplingTimes);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(lowerSamplingTimes, bt));
		String format = "\ngraph.infA<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphOnly(PrintWriter fout) {
		String decision = "FALSE";
		if (lexer.hasPlottingBlock()) decision = "TRUE";
		
		String format = "\ngraph.only<-%s\n";
		fout.write(String.format(format, decision));
	}
	
	private void writeGraphSpecification(PrintWriter fout) {
		writeSectionHeader(fout, "GRAPH SPECIFICATION OPTION");
		writeGraphLogical(fout);
		writeGraphLogicalSensitivity(fout);
		writeXAxesNames(fout);
		writeYAxesNames(fout);
		writeLogLogical(fout);
		writeGraphLowerSamplingTimes(fout);
		writeGraphUpperSamplingTimes(fout);
		writeYAxesRange(fout);
	}
	
	private void writeGraphUpperSamplingTimes(PrintWriter fout) {
		if (upperSamplingTimes == null) return;

		BinaryTree bt = tm.newInstance(upperSamplingTimes);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(upperSamplingTimes, bt));
		String format = "\ngraph.supA<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeIdenticalDose(PrintWriter fout) {
		BooleanValue value = null;
		if (usingIdenticalDose) value = new TrueBoolean();
		else value = new FalseBoolean();
		
		String flag = getSymbol(value);
		
		String format = "\ndose.identical<-%s %s Identical dose in each elementary design\n";
		fout.write(String.format(format, flag, comment_char));
	}
	
	private void writeIdenticalTimes(PrintWriter fout) {
		String value = parse(identicalTimes, tm.newInstance(identicalTimes)).trim();
		
		String format = "\nidentical.times<-%s %s Identical sampling times for each response.\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeInitialConditionsIdentical(PrintWriter fout) {
		BinaryTree bt = tm.newInstance(identicalInitialConditionsInEachDesign);
		
		String rValue = parse(identicalInitialConditionsInEachDesign, bt).trim();
		
		String format = "\ncondinit.identical<-%s\n";
		fout.write(String.format(format, rValue));
	}
	
	private void writeInitialConditionTime(PrintWriter fout) {
		if (initialConditionTime == null) return;
		
		BinaryTree bt = tm.newInstance(initialConditionTime);
		
		String value = parse(initialConditionTime, bt).trim();
		String format = "\ntime.condinit<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeLogLogical(PrintWriter fout) {
		String value = null;
		if (standardGraphic) value = "F";
		
		if (value == null) return;
		
		String format = "\nlog.logical<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeLowerA(PrintWriter fout) {
		if (lowerA == null) return;

		BinaryTree bt = tm.newInstance(lowerA);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(lowerA, bt));
		String format = "\nlowerA<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeLowerB(PrintWriter fout) {
		if (lowerB == null) return;

		BinaryTree bt = tm.newInstance(lowerB);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(lowerB, bt));
		String format = "\nlowerB<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeMajorOptions(PrintWriter fout) {
		writeSectionHeader(fout, "INPUT FILE FOR PFIM v" + lexer.getTarget().getVersion().getMajor());
		writeProjectName(fout);
		writeModelFile(fout);
		writeOutputFilename(fout);
		writeOutputFIMFilename(fout);
		writeFIM(fout);
		writePreviousFIM(fout);
		writeRun(fout);
		writeGraphOnly(fout);
		writeFimTypeOption(fout);
		writeNumberOfResponses(fout);
	}
	
	private void writeMaximumIterations(PrintWriter fout) {
		if (maximumIterations == null) return;

		String format = "\nmax.iter<-%s %s Parameter for initial simplex building\n";
		fout.write(String.format(format, maximumIterations, comment_char));
	}
	
	private void writeModelFile(PrintWriter fout) {
		String model_filename = getModelFunctionName(null) + "." + script_file_suffix;
		String format = "\nfile.model<-\"%s\" %s Model Filename\n";
		fout.write(String.format(format, model_filename, comment_char));
	}
	
	private void writeModelForm(PrintWriter fout) {
		String form = "AF";
		
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb.isODE()) form = "DE";
		
		String format = "\nmodelform<-\"%s\"\n";
		fout.write(String.format(format, form));
	}
	
	@Override
	protected void writeModelFunction(PrintWriter fout, StructuralBlock sb) throws IOException {
		if (fout == null) throw new NullPointerException();
		
		writeScriptHeader(fout, lexer.getModelFilename());
		if (!sb.isODE()) writeAnalyticalModelFunction(fout, sb);
		
		writePFIMProjectFile();
	}
	
	private void writeModelOption(PrintWriter fout) {
		writeSectionHeader(fout, "MODEL OPTION");
		writeModelForm(fout);
	}
	
	private void writeNNI(PrintWriter fout) {
		String value = parse(computePower, tm.newInstance(computeNNI)).trim();
		
		String format = "\ncompute.nni<-%s %s Compute No. individuals for given power of comparison test.\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeNumberOfOccassions(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "\nn_occ<-%s %s Number of Occassions\n";
		fout.write(String.format(format, this.numberOfOccassions, comment_char));
	}
	
	private void writeNumberOfResponses(PrintWriter fout) {
		String format = "\nnr<-%s %s Number of responses\n";
		
		fout.write(String.format(format, numberOfResponses, comment_char));
	}
	
	private void writeNumericalDerivative(PrintWriter fout) {
		// PharmML always has a structural model for this setting always 'True'.
		fout.write("\nNUM<-T\n");
	}
	
	private void writeOccassionParameters(PrintWriter fout) {
		StringBuffer outer_stmt = new StringBuffer();
		int i = 0;
		String format = "\"%s\"";
		for (String name : occ_parameter_map.keySet()) {
			if (i > 0) outer_stmt.append(",");
			
			ParametersAssociatedWithOccassion pawo = occ_parameter_map.get(name);
			int j = 0;
			StringBuffer inner_stmt = new StringBuffer(name);
			inner_stmt.append("=c(");
			for (CommonParameter p : pawo.params) {
				if (j > 0) inner_stmt.append(","); 
				inner_stmt.append(String.format(format, z.get(p.getSymbId())));
				j++;
			}
			inner_stmt.append(")");
			
			outer_stmt.append(inner_stmt);
			i++;
		}
		
		format = "\nparameter_occ.associated<-list ( %s )\n";
		fout.write(String.format(format, outer_stmt));
	}
	
	private void writeOmega(PrintWriter fout) {
		if (fout == null) return;
		if (omega == null) return;
		
		String value = parse(omega, tm.newInstance(omega)).trim();
		
		String format = "\nomega<-%s %s Inter-subject random effects\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeOnlyForOptimisation(PrintWriter fout) {
		writeSectionHeader(fout, "Optimisation Options");
		writeIdenticalTimes(fout);
		writeAlgorithmOption(fout);
		
		if(algo.equals(SIMPLEX)) writeSimplexSpecification(fout);
	}
	
	private void writeOutputFilename(PrintWriter fout) {
		String outputs_filename = pfimStdoutFilename + "." + script_file_suffix;
		String format = "\noutput<-\"%s\" %s file containing the PK or PD model\n";
		fout.write(String.format(format, outputs_filename, comment_char));
	}
	
	private void writeOutputFIMFilename(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "\noutputFIM<-\"%s\"\n";
		fout.write(String.format(format, outputFIMFilename));
	}
	
	private void writeParameterAssociated(PrintWriter fout) {
		CovariateBlock cb = lexer.getCovariateBlock();
		if (cb == null) return;
		
		Collection<CovariateParameterRef> refs = cb.getCovariateParameterRefs();
		if (refs == null) return;
		if (refs.size() == 0) return;
		
		List<String> decls = new ArrayList<String>();
		for (CovariateParameterRef ref : refs) decls.add(parse(ref, tm.newInstance(ref)).trim());
		
		StringBuffer stmt = new StringBuffer();
		int i = 0;
		for (String decl : decls) {
			if (i > 0) stmt.append(",");
			stmt.append(decl);
			i++;
		}
		
		String format = "\nparameter.associated<-list(%s)\n";
		fout.write(String.format(format, stmt));
	}
	
	private void writePFIMProjectFile() throws IOException {
		if (wrotePFIM_R) return;
		
		String cwd = lexer.getOutputDirectory();
		String outputFilepath = getPFIMProjectFilepath();
		
		PrintWriter fout = new PrintWriter(outputFilepath);
		
		writeScriptHeader(fout, null);
		
		// Filter the project template for project specific settings.
		for (int i = 0; i < pfimProjectTemplate.size(); i++) {
			String line = pfimProjectTemplate.get(i);
			if (line == null) continue;
			if (line.contains(PFIM_PROGRAM_DIR)) line = line.replace(PFIM_PROGRAM_DIR, programDirectory);
			else if (line.contains(CURRENT_WORKING_DIR)) line = line.replace(CURRENT_WORKING_DIR, cwd);
			
			fout.write(line + "\n");
		}
		fout.close();
		
		wrotePFIM_R = true;
	}
	
	private void writePowerAndNumberOfSubjects(PrintWriter fout) {
		writeSectionHeader(fout, "Power and number of subjects");
		writeAlpha(fout);
		writeComputePower(fout);
		writeNNI(fout);
		writeEquivalenceInterval(fout);
		writeComputePowerEquivalence(fout);
		writeComputeNNIEquivalence(fout);
		writeGivenPower(fout);
	}
	
	private void writePreviousFIM(PrintWriter fout) {
		String format = "\nprevious.FIM<-\"%s\"\n";
		fout.write(String.format(format, previousFIM));
	}
	
	private void writePrintIterations(PrintWriter fout) {
		if (printIterations == null) return;

		String value = parse(printIterations, tm.newInstance(printIterations)).trim();
		String format = "\niter.print<-%s %s Print Iterations\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeProjectName(PrintWriter fout) {
		String format = "\nproject<-\"%s\" %s Project Name\n";
		fout.write(String.format(format, lexer.getModelName(), comment_char));
	}
	
	private void writeProportions(PrintWriter fout) {
		String option = null;
		if (usingSubjects) option = "1"; // No. Subject
		else option = "2"; // Using proportions
		
		String format = "\nsubjects.input<-%s %s Subjects Input\n";
		fout.write(String.format(format, option, comment_char));
	}
	
	private void writeProtA(PrintWriter fout) {
		if (protA.isEmpty()) return;
		
		List<String> stmts = new ArrayList<String>();  
		for (Object samplingTimes : protA) {
			if (samplingTimes == null) continue;
			BinaryTree bt = tm.newInstance(samplingTimes);
			lexer.updateNestedTrees();
			stmts.add(parse(protA, bt).trim());
		}
		if (stmts.isEmpty()) return;
		
		StringBuffer list = new StringBuffer("list(");
		int i = 0;
		for (String clause: stmts) {
			if (i > 0) list.append(",");
			list.append(clause);
			i++;
		}
		list.append(")");
		
		String format = "\nprotA<-%s %s Sampling Times\n";
		fout.write(String.format(format, list, comment_char));
	}
	
	@Override
	protected void writePurgeCommand(PrintWriter fout) {}
	
	private void writeRandomEffectModel(PrintWriter fout) {
		String format = "\nTrand<-%s; %s Random Effect model\n";
		fout.write(String.format(format, Trand, comment_char));
	}
	
	private void writeRelativeConverganceTolerance(PrintWriter fout) {
		if (maximumIterations == null) return;

		String format = "\nRctol<-%s %s Relative convergence tolerance\n";
		fout.write(String.format(format, relativeConvergenceTolerance, comment_char));
	}
	
	private void writeResidualError(PrintWriter fout) {
		String format = "\n%s Standard deviation of residual error\n";
		
		fout.write(String.format(format, comment_char));
		
		format = "%s<-%s\n";
		fout.write(String.format(format, "sig.interA", sigInterA));
		fout.write(String.format(format, "sig.slopeA", sigSlopeA));
	}
	
	private void writeRun(PrintWriter fout) {
		if (fout == null) return;
		
		String run = null;
		
		if (lexer.hasSimulation()) run = "EVAL";
		else if (lexer.hasEstimation()) run = "OPT";
		
		if (run == null) throw new IllegalStateException("Unable to determine the 'RUN' option.");
		
		String format = "\nrun<-\"%s\"\n";
		fout.write(String.format(format, run));
	}
	
	@Override
	protected void writeScriptHeader(PrintWriter fout, String model_file) throws IOException {
		if (fout == null) return;

		String format = "%s Script generated by the PFIM Converter\n";
		fout.write(String.format(format, comment_char));

		format = "%s Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM (2016)\n";
		fout.write(String.format(format, comment_char));
		
		format = "%s Converter Version: %s\n";
		fout.write(String.format(format, comment_char, lexer.getConverterVersion()));
		
		format = "%s Source: %s\n";
		fout.write(String.format(format, comment_char, lexer.getSource()));
		
		format = "%s Target: %s\n";
		fout.write(String.format(format, comment_char, lexer.getTarget()));

		format = "%s Run ID: %s\n";
		fout.write(String.format(format, comment_char, run_id));
		
		String title = lexer.getModelName();
		if (title != null) {
			format = "%s Model: %s\n";
			fout.write(String.format(format, comment_char, title));
		}

		if (model_file != null) {
			format = "%s File: %s\n";
			fout.write(String.format(format, comment_char, model_file));
		}

		format = "%s Dated: %s\n\n";
		fout.write(String.format(format, comment_char, new Date()));
	}
	
	@Override
	protected void writeScriptLibraryReferences(PrintWriter fout) throws IOException {
		if (use_compiled_structural_model) fout.write("library(compiler)\n");
		
		String format = "\nsource('%s')\n";
		fout.write(String.format(format, this.getPFIMProjectFilepath()));
		fout.write("PFIM()\n");
	}
	
	private void writeSectionHeader(PrintWriter fout, String headerTitle) {
		String divider = "\n%s===============================================\n";
		String format = "%s\t%s";
		
		fout.write(String.format(divider, comment_char));
		fout.write(String.format(format, comment_char, headerTitle));
		fout.write(String.format(divider, comment_char));
	}
	
	private void writeSimplexParameter(PrintWriter fout) {
		if (simplexParameter == null) return;

		String format = "\nsimplex.parameter<-%s %s Parameter for initial simplex building\n";
		fout.write(String.format(format, simplexParameter, comment_char));
	}
	
	private void writeSimplexSpecification(PrintWriter fout) {
		writeSectionHeader(fout, "SIMPLEX SPECIFICATION");
		writeSubjectsOptimisation(fout);
		writeLowerA(fout);
		writeUpperA(fout);
		writeLowerB(fout);
		writeUpperB(fout);
		writeDeltaTime(fout);
		writePrintIterations(fout);
		writeSimplexParameter(fout);
		writeMaximumIterations(fout);
		writeRelativeConverganceTolerance(fout);
	}
	
	private void writeSolverSettings(PrintWriter fout) {
		String format = "\n%s Error tolerance for solving differential equations\n";
		fout.write(String.format(format, comment_char));
		
		format = "%s<-%s\n";
		fout.write(String.format(format, "RtolEQ", rtol));
		fout.write(String.format(format, "AtolEQ", atol));
		fout.write(String.format(format, "Hmax", hmax));
	}
	
	private void writeStdin() throws IOException {
		if (writtenSTDIN) return;
		
		String outFilepath = getStdinFilepath();
		
		readOperations();

		tm = lexer.getTreeMaker();
		PrintWriter fout = new PrintWriter(outFilepath);
		
		writeScriptHeader(fout, lexer.getModelFilename());
		writeMajorOptions(fout);
		writeModelOption(fout);
		writeAnalyticalModelOption(fout);
		writeDifferentialEquationOption(fout);
		writeEffectsAndErrorModel(fout);
		writeCovariateModel(fout);
		writeCovariatesChangingWithOccassion(fout);
		writePowerAndNumberOfSubjects(fout);
		writeOnlyForOptimisation(fout);
		writeGraphSpecification(fout);
		
		fout.close();
		
		writtenSTDIN = true;
	}
	
	private void writeSubjects(PrintWriter fout) {
		if (!usingSubjects) return;
		if (numberOfSubjects == 0) return;
		
		String format = "\nsubjects<-c(%s) %s No. Subjects\n";
		fout.write(String.format(format, numberOfSubjects, comment_char));
	}
	
	private void writeSubjectsOptimisation(PrintWriter fout) {
		String value = parse(optimisationOfProportionsOfSubjects, tm.newInstance(optimisationOfProportionsOfSubjects)).trim();
		String format = "\nsubjects.opt<-%s %s OPT on proportions of subjects\n";
		fout.write(String.format(format, value, comment_char));
	}
	
	private void writeUpperA(PrintWriter fout) {
		if (upperA == null) return;

		BinaryTree bt = tm.newInstance(upperA);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(upperA, bt));
		String format = "\nupperA<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeUpperB(PrintWriter fout) {
		if (upperB == null) return;

		BinaryTree bt = tm.newInstance(upperB);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(upperB, bt));
		String format = "\nupperB<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeUsingCovariates(PrintWriter fout) {
		String option = "FALSE";
		if (usingCovariateModel) option = "TRUE";
		
		String format = "\ncovariate.model<-%s %s Add covariate to the model\n";
		fout.write(String.format(format, option, comment_char));
	}
	
	private void writeXAxesNames(PrintWriter fout) {
		if (xAxesNames == null) return;

		BinaryTree bt = tm.newInstance(xAxesNames);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(xAxesNames, bt));
		String format = "\nnames.datax<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeYAxesNames(PrintWriter fout) {
		if (yAxesNames == null) return;

		BinaryTree bt = tm.newInstance(yAxesNames);
		lexer.updateNestedTrees();
		
		String value = stripOuterBrackets(parse(yAxesNames, bt));
		String format = "\nnames.datay<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeYAxesRange(PrintWriter fout) {
		String value = "NULL";
		
		if (yAxesRange != null) {
			BinaryTree bt = tm.newInstance(yAxesRange);
			lexer.updateNestedTrees();
			value = stripOuterBrackets(parse(yAxesRange, bt));
		}
		
		String format = "\ny.rangeA<-%s\n";
		fout.write(String.format(format, value));
	}
}

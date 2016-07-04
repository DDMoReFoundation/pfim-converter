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

import static crx.converter.engine.PharmMLTypeChecker.isBinaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isConstant;
import static crx.converter.engine.PharmMLTypeChecker.isDerivative;
import static crx.converter.engine.PharmMLTypeChecker.isFalse;
import static crx.converter.engine.PharmMLTypeChecker.isFunctionCall;
import static crx.converter.engine.PharmMLTypeChecker.isIndependentVariable;
import static crx.converter.engine.PharmMLTypeChecker.isInt;
import static crx.converter.engine.PharmMLTypeChecker.isJAXBElement;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isLogicalBinaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isLogicalUnaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isParameterEstimate;
import static crx.converter.engine.PharmMLTypeChecker.isPiecewise;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;
import static crx.converter.engine.PharmMLTypeChecker.isRandomVariable;
import static crx.converter.engine.PharmMLTypeChecker.isReal;
import static crx.converter.engine.PharmMLTypeChecker.isRhs;
import static crx.converter.engine.PharmMLTypeChecker.isSequence;
import static crx.converter.engine.PharmMLTypeChecker.isString;
import static crx.converter.engine.PharmMLTypeChecker.isStructuredError;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isTrue;
import static crx.converter.engine.PharmMLTypeChecker.isUnaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isVariableReference;
import static crx.converter.engine.PharmMLTypeChecker.isVector;
import static inserm.converters.pfim.OptimisationAlgorithm.FEDOROV_WYNN;
import static inserm.converters.pfim.OptimisationAlgorithm.SIMPLEX;
import static inserm.converters.pfim.SettingLabel.GRAPH_SUPA;
import static inserm.converters.pfim.SettingLabel.OUTPUT;
import static inserm.converters.pfim.SettingLabel.OUTPUT_FIM;
import static inserm.converters.pfim.SettingLabel.PROJECT;
import inserm.converters.pfim.parts.OptimalDesignStepImpl;
import inserm.converters.pfim.parts.ParameterBlockImpl;
import inserm.converters.pfim.parts.TrialDesignBlockImpl;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.Accessor;
import crx.converter.engine.ConversionDetail_;
import crx.converter.engine.SymbolReader.ModifiedSymbol;
import crx.converter.engine.common.BaseParser;
import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.engine.common.InterventionSequenceRef;
import crx.converter.engine.common.SimulationOutput;
import crx.converter.spi.OptimalDesignLexer;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.blocks.TrialDesignBlock2;
import crx.converter.spi.blocks.VariabilityBlock;
import crx.converter.spi.steps.OptimalDesignStep_;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.Node;
import crx.converter.tree.TreeMaker;
import eu.ddmore.convertertoolbox.api.response.ConversionDetail;
import eu.ddmore.libpharmml.dom.IndependentVariable;
import eu.ddmore.libpharmml.dom.commontypes.BooleanValue;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.InitialCondition;
import eu.ddmore.libpharmml.dom.commontypes.IntValue;
import eu.ddmore.libpharmml.dom.commontypes.LevelReference;
import eu.ddmore.libpharmml.dom.commontypes.OidRef;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.Sequence;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.TrueBoolean;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.commontypes.Vector;
import eu.ddmore.libpharmml.dom.commontypes.VectorElements;
import eu.ddmore.libpharmml.dom.commontypes.VectorValue;
import eu.ddmore.libpharmml.dom.maths.Constant;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType;
import eu.ddmore.libpharmml.dom.maths.FunctionCallType.FunctionArgument;
import eu.ddmore.libpharmml.dom.maths.LogicBinOp;
import eu.ddmore.libpharmml.dom.maths.LogicUniOp;
import eu.ddmore.libpharmml.dom.maths.Piece;
import eu.ddmore.libpharmml.dom.maths.Piecewise;
import eu.ddmore.libpharmml.dom.modeldefn.Distribution;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.ObservationError;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomEffect;
import eu.ddmore.libpharmml.dom.modeldefn.ParameterRandomVariable;
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel.LinearCovariate;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredModel.LinearCovariate.PopulationValue;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError.ErrorModel;
import eu.ddmore.libpharmml.dom.modellingsteps.FIMtype;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;
import eu.ddmore.libpharmml.dom.probonto.DistributionParameter;
import eu.ddmore.libpharmml.dom.probonto.ParameterName;
import eu.ddmore.libpharmml.dom.probonto.ProbOnto;
import eu.ddmore.libpharmml.dom.trialdesign.ArmDefinition;
import eu.ddmore.libpharmml.dom.trialdesign.Observation;
import eu.ddmore.libpharmml.dom.trialdesign.SingleDesignSpace;
import eu.ddmore.libpharmml.dom.uncertml.VarRefType;

/**
 * PFIM R-based code generator.
 */
public class Parser extends BaseParser {
	private final static String directory_program_replacement_label = "DIRECTORY_PROGRAM_REPLACE";
	private final static String directory_replacement_label = "DIRECTORY_REPLACE";
	private static final String MODEL_FILESTEM = "model";
	private static final String pfimFIMFilename = "FIM.txt";
	private static final String pfimProjectFilename = "PFIM";
	private static final String pfimStdinFilename = "stdin";
	private static final String pfimStdoutFilename = "stdout.out";
	private static final String PREFERRED_SEPERATOR = "/";
	
	private static String toAlphabetic(int i) {
	    if(i < 0) return "-"+ toAlphabetic(-i-1);
	    
	    int quot = i / 26;
	    int rem = i % 26;
	    char letter = (char)((int) 'A' + rem);
	    if (quot == 0) return ""+letter;
	    else return toAlphabetic(quot-1) + letter;
	}
	
	private OptimisationAlgorithm algorithm = FEDOROV_WYNN;
	private double alpha = 0.05;
	private String analytical_model_dose_variable = "X";
	private Double atol = 1E-08;
	private boolean computeNNI = false;
	private boolean computeNNIEquivalence = false;
	private boolean computePower = false;
	private boolean computePowerEquivalence = false;
	private Object ctx = new Object();
	private double current_offset = 0.0;
	private double deltaTime = 0.0;
	private Map<ArmDefinition, String> dosing_vector_map = new HashMap<ArmDefinition, String>();
	private FIMtype fim_type = FIMtype.P;
	private double givenPower = 0.9;
	private boolean identicalTimes = true;
	private String leftArrayBracket = null;
	private int maximumIterations = 5000;
	private Map<Observation, List<Double>> ob_sampling_map = new HashMap<Observation, List<Double>>();
	private OutputFIMFormat option = OutputFIMFormat.BLOCK_DIAGONAL_FIM;
	private boolean optProportionsOfsubjects = false;
	private String output_state_vector_symbol = "yd";
	private String outputFIMFilename = pfimFIMFilename;
	private String param_model_symbol  = null;
	private List<String> pfimProjectTemplate = new ArrayList<String>();
	private String previousFIM = "";
	private boolean printIterations = true;
	private String programDirectory = ".";
	private Properties props = null;
	private boolean record_vector_values = false;
	private List<Double> recorded_vector_values = new ArrayList<Double>();
	private Map<IndividualParameter, PopulationParameter> referenced_params_map = new HashMap<IndividualParameter, PopulationParameter>();
	private double relativeConvergenceTolerance = 1E-6;
	private String rightArrayBracket = null;
	private Double rtol = 1E-08;
	private List<Double> secondary_recorded_values = null;
	private int simplexParameter = 20;
	private SettingReader sr = null;
	private String state_vector_symbol = null;
	private String stdoutFilename = pfimStdoutFilename;
	private RandomEffectModelOption trand = RandomEffectModelOption.EXPONENTIAL;
	
	/**
	 * Constructor
	 * @throws IOException
	 */
	public Parser() throws IOException {
		super();
		useConversativeBinaryOperatorFormats();
		loadPFIMTemplate();
		
		props = new Properties();
		props.load(getClass().getResourceAsStream("Parser.properties"));
		
		comment_char = "#";
		script_file_suffix = "R";
		objective_dataset_file_suffix = "csv";
		output_file_suffix = "csv";
		solver = "ode";
		leftArrayBracket = "[";
		rightArrayBracket = "]";
		
		// PFIM Specific symbols and settings.
		param_model_symbol = "p";
		state_vector_symbol = "y";
		programDirectory = props.getProperty("pfimProgramDirectory");
		
		setReferenceClass(getClass());
		init();
		
		// Read directory path via the shell, ignoring the relative path
		// specified in the config file.
		if (System.getProperty("PFIM_PROG_DIR") != null) programDirectory = System.getProperty("PFIM_PROG_DIR");
	}
	
	private String cat(List<String> values) {
		if (values == null) return "c()";
		StringBuffer s = new StringBuffer("c(");
		int i = 0;
		for (String value : values) {
			if (value == null) continue;
			if (i > 0) s.append(",");
			s.append(value);
			i++;
		}
		s.append(")");
		return s.toString();
	}
	
	private String cat_(List<Double> values) {
		if (values == null) return "c()";
		StringBuffer s = new StringBuffer("c(");
		int i = 0;
		for (Double value : values) {
			if (value == null) continue;
			if (i > 0) s.append(",");
			s.append(value);
			i++;
		}
		s.append(")");
		return s.toString();
	}
	
	void createSettingReader() {
		Converter c = (Converter) lexer;
		
		sr = new SettingReader();
		sr.setLexer(c);
		sr.setParser(c.getParser());
		sr.setStep(c.getOptimalDesignStep().getStep());
		sr.readSettings();
	}
	
	private String doArrayAccess(String variableName, Integer idx) {
		String format = "%s%s%s%s";
		
		if (variableName == null || idx == null) throw new NullPointerException("Required array access data is NULL");
		
		return String.format(format, variableName, leftArrayBracket, idx, rightArrayBracket);
	}
	
	private String doBigInteger(BigInteger i) { return i.toString(); }
	
	private String doConstant(Constant c) {
		String symbol = unassigned_symbol;
		
		String op = c.getOp();
		if (op.equalsIgnoreCase("notanumber")) symbol = "NaN";
		else if (op.equalsIgnoreCase("pi")) symbol = "pi";
		else if (op.equalsIgnoreCase("exponentiale")) symbol = "exp(1)";
		else if (op.equalsIgnoreCase("infinity")) symbol = "Inf";
	
		return symbol;
	}
	
	private String doDerivative(DerivativeVariable s) {
		Integer idx = lexer.getStateVariableIndex(s.getSymbId());
		String format = "%s%s";
		return String.format(format, output_state_vector_symbol, idx);
	}
	
	private String doElement(JAXBElement<?> element) {
		String symbol = unassigned_symbol;
		
		Object value = element.getValue();
		if (value != null) {
			if (isBinaryOperation(value) || isUnaryOperation(value)) 
				symbol = parseRawEquation(value);
			else 
				symbol = getSymbol(element.getValue());
		}
		
		return symbol;
	}
	
	private String doFalse() { return "FALSE"; }
	
	private String doFunctionCall(FunctionCallType call) {
		if (call == null) throw new NullPointerException("A function call definition is null.");
		
		ArrayList<String> arg_symbols = new ArrayList<String>();
		Object context = new Object();
		for (FunctionArgument arg : call.getListOfFunctionArgument()) {
			if (arg == null) continue;
			
			String arg_symbol = "0.0";
			TreeMaker tm = lexer.getTreeMaker();
			if(isPiecewiseFunctionArgument(arg)) arg_symbol = doPiecewiseForFunctionArgument(arg);
			else {
				BinaryTree bt = null;
				if (!lexer.hasStatement(arg)) {
					bt = tm.newInstance(arg);
					lexer.updateNestedTrees();
				} else
					bt = lexer.getStatement(arg);
				
				arg_symbol = parse(context, bt).trim();	
			}
			
			arg_symbols.add(arg_symbol);
		}

		StringBuffer args_string = new StringBuffer();
		int i = 0;
		for (String symbol : arg_symbols) {
			if (i > 0)
				args_string.append(", ");
			args_string.append(symbol);
			i++;
		}
		
		String format = "%s(%s)";
		return String.format(format, z.get(call), args_string);
	}
	
	private String doIndependentVariable(IndependentVariable v) {
		return z.get(v.getSymbId());
	}
	
	private String doIndividualParameterAssignment(IndividualParameterAssignment ipa) {
		ParameterBlockImpl pb = (ParameterBlockImpl) lexer.getParameterBlock();
		return doArrayAccess(param_model_symbol, pb.getParameterIndex(ipa.parameter.getSymbId()));
	}
	
	private String doInt(IntValue i) {
		return i.getValue().toString();
	}
	
	private String doJavaBoolean(Boolean value) {
		if (value) return doTrue();
		else return doFalse();
	}
	
	/**
	 * Write code for local variable reference.
	 * @param v Local Variable
	 * @return String
	 */
	protected String doLocalVariable(VariableDefinition v) {
		return z.get(v.getSymbId());
	}
	
	protected String doLogicalBinaryOperator(LogicBinOp l_b_op) {
		return getLogicalOperator(l_b_op.getOp());
	}
	
	protected String doLogicalUnaryOperator(LogicUniOp u_b_op) {
		return getLogicalOperator(u_b_op.getOp());
	}
	
	/**
	 * Generate code for a parameter reference.
	 * @param p Parameter
	 * @return String
	 */
	protected String doParameter(PopulationParameter p) { return z.get(p); }
	
	protected String doParameterEstimate(ParameterEstimate pe) {
		return unassigned_symbol;
	}
	
	public String doPiecewise(Piecewise pw) { return unassigned_symbol; }
	
	protected String doPiecewiseForFunctionArgument(FunctionArgument arg) {
		if (!isPiecewiseFunctionArgument(arg)) return unassigned_symbol;
		
		Piecewise pw = arg.getAssign().getPiecewise();
		List<Piece> pieces = pw.getListOfPiece();
		if (pieces.isEmpty()) throw new IllegalStateException("The function argument has no expected piecewise statement.");
		
		TreeMaker tm = lexer.getTreeMaker();
		String format = "(%s)";
		return String.format(format, parse(new Object(), tm.newInstance(pieces.get(0).getCondition())));
	}
	
	private String doReal(RealValue r) { return Double.toString(r.getValue()); }
	
	private String doRhs(Rhs eq) {
		TreeMaker tm = lexer.getTreeMaker();
		return parse(new Object(), tm.newInstance(eq));
	}
	
	private String doSequence(Sequence o) {
		String symbol = unassigned_symbol;
		
		Sequence seq = (Sequence) o;
		Rhs begin = seq.getBegin();
		Rhs end = seq.getEnd();
		Rhs repetitions = seq.getRepetitions();
		Rhs stepSize = seq.getStepSize();
		
		BinaryTree btBegin = null, btEnd = null, btRepetitions = null, btStepSize = null;
		if (begin == null) throw new IllegalStateException("The required Sequence.begin field is not assigned.");
		
		if (begin != null) btBegin = lexer.getStatement(begin);
		if (end != null) btEnd = lexer.getStatement(end);
		if (repetitions != null) btRepetitions = lexer.getStatement(repetitions);
		if (stepSize != null) btStepSize = lexer.getStatement(stepSize);
		
		String strBegin = null, strEnd = null, strRepetitions = null, strStepSize = null;
		
		if (btBegin != null) strBegin = parse(seq, btBegin);
		if (btEnd != null) strEnd = parse(seq, btEnd);	
		if (btRepetitions != null) strRepetitions = parse(seq, btRepetitions);
		if (btStepSize != null)	strStepSize = parse(seq, btStepSize);
			
		// Default value in case the conditional logic fails or the PharmML spec changes.
		symbol = "c(0 0)";
		if (strBegin != null && strEnd != null && strStepSize != null) {
			String format = "seq(%s, %s, by=%s)";
			symbol = String.format(format, strBegin, strEnd, strStepSize);
		} else if (strBegin != null && strEnd != null && strRepetitions != null) {
			String format = "seq(%s,%s,length=%s)";
			symbol = String.format(format, strBegin, strEnd, strRepetitions);
		} else if (strBegin != null && strStepSize != null && strRepetitions != null) {
			String format = "seq(%s, (%s*%s), by=%s)";
			symbol = String.format(format, strBegin, strStepSize, strRepetitions, strStepSize);
		} else if (strBegin != null && strEnd != null) {
			String format = "%s:%s";
			symbol = String.format(format, strBegin, strEnd);
		}
		
		return symbol;
	}
	
	private String doString(String v) { return v; }
	
	private String doStringValue(StringValue sv) {
		String format = "'%s'";
		return String.format(format, sv.getValue());
	}
	
	private String doSymbolRef(SymbolRef s) {
		String symbol = s.getSymbIdRef();
		
		Accessor a = lexer.getAccessor();
		PharmMLRootType element = a.fetchElement(s);
		
		if (lexer.isStateVariable(symbol) || isDerivative(element)) 
			symbol = doArrayAccess(state_vector_symbol, lexer.getStateVariableIndex(symbol));
		else {
			if (lexer.isRemoveIllegalCharacters()) {
				ModifiedSymbol result = z.removeIllegalCharacters(s, symbol);
				if (result.isModified()) symbol = result.modified_value;
			}
			
			if (lexer.isFilterReservedWords()) {
				if (z.isReservedWord(symbol)) {
					symbol = z.replacement4ReservedWord(s.getSymbIdRef());
					if (symbol == null) 
						throw new NullPointerException("Replacement symbol for reserved word (symbIdRef=('" + s.getSymbIdRef() + "') undefined.");
					
					ModifiedSymbol result = new ModifiedSymbol(element, s.getSymbIdRef(), symbol);
					if (result.isModified()) z.add(result);
				}
			}
		}
		
		return symbol;
	}
	
	private String doTrue() { return "TRUE"; }
	
	private String doVarRef(VarRefType ref) {
		String symbol = unassigned_symbol;
		
		PharmMLRootType element = lexer.getAccessor().fetchElement(ref);
				
    	if (element == null) throw new NullPointerException("Variable reference not found (" + ref.getVarId() + ")");
    	symbol = getSymbol(element);
		
		return symbol;
	}
	
	private String doVector(Vector v) {
		String symbol = unassigned_symbol;
		
		List<String> values = new ArrayList<String>();
		
		VectorElements elements = v.getVectorElements();
		if (elements == null) throw new NullPointerException("Vector elements are NULL.");
		
		if (elements.getListOfElements().size() == 1) {
			VectorValue value = elements.getListOfElements().get(0);
			if (value != null) return "0.0";
		}
		
		for (VectorValue value : elements.getListOfElements()) {
			if (value == null) continue;
			if (lexer.hasStatement(value)) {
				String stmt = parse(v, lexer.getStatement(value));
				stmt = stripOuterBrackets(stmt);
				values.add(stmt);
			}
		}
		
		if (record_vector_values) recordVectorValues(values);
		
		StringBuilder st = new StringBuilder();
		st.append("c(");
		int count = 0;
		for (String value : values) {
			if (count > 0) st.append(",");
			st.append(value.trim());
			count++;
		}
		st.append(")");
		
		symbol = st.toString();
		return stripOuterBrackets(symbol);
	}
	
	private String getBatchFilepath() {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + "run.bat";
	}
	
	private String getModelFilename() {
		return  MODEL_FILESTEM + "." + script_file_suffix;
	}
	
	@Override
	public String getModelFunctionFilename(String output_dir, StructuralBlock sb) {
		return String.format("%s%s%s.%s", output_dir, File.separator, MODEL_FILESTEM, script_file_suffix);
	}
	
	private String getParameterValueFromEstimate(SymbolRef ref) {
		String slope = "0.0";
		if (ref == null) return slope;
		
		Accessor a = lexer.getAccessor();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		
		PharmMLRootType element = a.fetchElement(ref);
		if (!isPopulationParameter(element)) throw new NullPointerException("The proportional error component is not a numeric constant.");
		PopulationParameter p = (PopulationParameter) element;
		ParameterEstimate pe = step.getParameterEstimate(p);
		if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
		slope = parse(ctx, lexer.getStatement(pe.getInitialEstimate())).trim();
		return slope;
	}
	
	private String getPFIMProjectFilepath() throws IOException {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + pfimProjectFilename + "." + script_file_suffix;
	} 
	
	private String getStdinFilepath() {
		String cwd = lexer.getOutputDirectory();
		String filename = pfimStdinFilename;
		return cwd + PREFERRED_SEPERATOR + filename + "." + script_file_suffix;
	} 
	
	@Override
	public String getSymbol(Object o) {
		String symbol = unassigned_symbol;
		
		if (isSymbolReference(o)) symbol = doSymbolRef((SymbolRef) o); 
		else if (isDerivative(o)) symbol = doDerivative((DerivativeVariable) o);
		else if (isPopulationParameter(o)) symbol = doParameter((PopulationParameter) o);
		else if (isLocalVariable(o)) symbol = doLocalVariable((VariableDefinition) o);
		else if (isString_(o))  symbol = doString((String) o);
		else if (isReal(o)) symbol = doReal((RealValue) o);
		else if (isFalse(o)) symbol = doFalse();
		else if (isTrue(o)) symbol = doTrue();
		else if (isString(o)) symbol = doStringValue((StringValue) o);
		else if (isInt(o)) symbol = doInt((IntValue) o);
		else if (isConstant(o)) symbol = doConstant((Constant) o);
		else if (isIndependentVariable(o)) symbol = doIndependentVariable((IndependentVariable) o);	
		else if (isFunctionCall(o)) symbol = doFunctionCall((FunctionCallType) o);
		else if (isSequence(o)) symbol = doSequence((Sequence) o); 
		else if (isVector(o)) symbol = doVector((Vector) o);
		else if (isPiecewise(o)) symbol =  doPiecewise((Piecewise) o);
		else if (isLogicalBinaryOperation(o)) symbol = doLogicalBinaryOperator((LogicBinOp) o);
	    else if (isLogicalUnaryOperation(o)) symbol = doLogicalUnaryOperator((LogicUniOp) o);
	    else if (isParameterEstimate(o)) symbol = doParameterEstimate((ParameterEstimate) o);
	    else if (isJAXBElement(o)) symbol = doElement((JAXBElement<?>) o);	
	    else if (isVariableReference(o)) symbol = doVarRef((VarRefType) o);	
	    else if (o instanceof IndividualParameterAssignment) symbol = doIndividualParameterAssignment((IndividualParameterAssignment) o);	
	    else if (isBigInteger(o)) symbol = doBigInteger((BigInteger) o);
	    else if (isRhs(o)) symbol = doRhs((Rhs) o);
	    else if (o instanceof Boolean) symbol = doJavaBoolean((Boolean) o);
	    else 
	    {
			String format = "WARNING: Unknown symbol, %s\n";
			String value = "null";
			if (o != null) value = o.toString();
			String msg = String.format(format, value);
			ConversionDetail detail = new ConversionDetail_();
			detail.setSeverity(ConversionDetail.Severity.WARNING);
			detail.addInfo("warning", msg);
			
			System.err.println(msg); 
		}
		
		return symbol;
	}
	
	public void initialise() {
		// Check if replicate setting via PharmML tags.
		if (!lexer.hasExternalDatasets()) lexer.checkForTrialDesignIfRandomisedModel();
		
		lexer.setRemoveIllegalCharacters(true);
		z.setReplacementCharacter('_');
		z.setIllegalCharacters(new char [] {'~', '@', '+', '*', '-', '/', '$', '!', '.', '(', ')', '[', ']', ',', '#', '%', "\n".charAt(0), " ".charAt(0)});
		
		lexer.setFilterReservedWords(true);
		try {
			z.loadReservedWords();
		} catch (Exception e) {
			e.printStackTrace(System.err);
			System.err.println("WARNING: Failed to read reserved word map for SymbolReader.");
		}
		
		lexer.setSaveRenamedSymbolList(true);
		setRunId("call_run");
	}
	
	private boolean is_parameter_scope(ParameterRandomVariable eta) {
		if (eta == null) return false;
		
		LevelReference lref = eta.getListOfVariabilityReference().get(0);
		if (lref != null) {
			VariabilityBlock vb = lexer.getVariabilityBlock(lref.getSymbRef());
			if (vb != null) return vb.isParameterVariability();
		}
		
		return false;
	}
	
	private boolean isPiecewiseFunctionArgument(FunctionArgument arg) {
		if (arg == null) return false;
		if (arg.getAssign() != null) {
			Rhs rhs = arg.getAssign();
			return isPiecewise(rhs.getContent());
		}
		
		return false;
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
	
	private void processErrorModel(int idx, ObservationBlock ob, PrintWriter fout) {
		if (ob == null || fout == null) return;

		ObservationError error_model = ob.getObservationError();
		if (!isStructuredError(error_model)) return;
		StructuredObsError serr = (StructuredObsError) ob.getObservationError();
		ErrorModel model = serr.getErrorModel();
		if (model == null) return;

		Rhs assign = model.getAssign();
		if (assign == null) return;

		Object content = assign.getContent();
		if (!isFunctionCall(content)) throw new IllegalStateException("An error model does not define the expected function call.");

		FunctionCallType call = (FunctionCallType) content;
		String func_name = call.getSymbRef().getSymbIdRef();
		if (!SupportedErrorModel.contains(func_name)) throw new UnsupportedOperationException("Error model not supported (type='" + func_name + "')");

		String inter = "0.0", slope = "0.0";		
		SupportedErrorModel model_flag = SupportedErrorModel.fromValue(func_name);
		
		final String proportional = "proportional", additive= "additive";
		if (SupportedErrorModel.PROPORTIONAL.equals(model_flag) || 
			SupportedErrorModel.PROPORTIONAL_ERROR.equals(model_flag)) {
			for (FunctionArgument arg : call.getListOfFunctionArgument()) {
				if (proportional.equals(arg.getSymbId())) {
					if (arg.getSymbRef() != null) slope = getParameterValueFromEstimate(arg.getSymbRef());
					else if (arg.getAssign() != null) {
						content = arg.getAssign().getContent();
						if (isSymbolReference(content)) slope = getParameterValueFromEstimate((SymbolRef) content);
					}
					
					break;
				}
			}
		} else if (SupportedErrorModel.COMBINED_ERROR.equals(model_flag)) {
			for (FunctionArgument arg : call.getListOfFunctionArgument()) {
				if (proportional.equals(arg.getSymbId())) {
					if (arg.getSymbRef() != null) slope = getParameterValueFromEstimate(arg.getSymbRef());
					else if (arg.getAssign() != null) {
						content = arg.getAssign().getContent();
						if (isSymbolReference(content)) slope = getParameterValueFromEstimate((SymbolRef) content);
					}
				} else if (additive.equals(arg.getSymbId())) {
					if (arg.getSymbRef() != null) inter = getParameterValueFromEstimate(arg.getSymbRef());
					else if (arg.getAssign() != null) {
						content = arg.getAssign().getContent();
						if (isSymbolReference(content)) inter = getParameterValueFromEstimate((SymbolRef) content);
					}
				}
			}
		}
		
		String label = toAlphabetic(idx);
		String format = "sig.inter%s<-%s\n";
		fout.write(String.format(format, label, inter));
		
		format = "sig.slope%s<-%s\n";
		fout.write(String.format(format, label, slope));
	}
	
	// Assuming a single random effect per IP instance.
	private String readOmega(ParameterRandomEffect rve) {
		String omega = "0.0";
		if (rve == null) return omega;
		
		Accessor a = lexer.getAccessor();
		SymbolRef ref = rve.getSymbRef().get(0);
		if (ref == null) return omega;
		
		
		PharmMLRootType element = a.fetchElement(ref);
		if (!isRandomVariable(element)) return omega;
		ParameterRandomVariable rv = (ParameterRandomVariable) element;
		if (!is_parameter_scope(rv)) 
			throw new IllegalStateException("Random variable does not have the required parameter variability scope (name='" + rv.getSymbId() + "')");
		
		// Just assuming probonto and normal usage, nothing else.
		Distribution dist = rv.getDistribution();
		if (dist == null) return omega;
		
		ProbOnto probonto = dist.getProbOnto();
		if (probonto == null) return omega;
		
		for (DistributionParameter dp : probonto.getListOfParameter()) {
			if (dp == null) continue;
			ParameterName name = dp.getName();
			if (name == null) continue;
			if (ParameterName.STDEV.equals(name)) {
				Object content = dp.getAssign().getContent();
				if (!isSymbolReference(content)) return omega;
				
				ref = (SymbolRef) content;
				element = a.fetchElement(ref);
				if (!isPopulationParameter(element)) return omega;
				PopulationParameter p = (PopulationParameter) element;
				
				OptimalDesignLexer c = (OptimalDesignLexer) lexer;
				OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
				ParameterEstimate pe = step.getParameterEstimate(p);
				if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
				omega = parse(ctx, lexer.getStatement(pe.getInitialEstimate())).trim();
			}
		}
		
		return omega;
	}

	private void recordVectorValues(List<String> values) {
		if (values == null) return;
		for (String value : values) {
			if (value == null) continue;
			try {
				Double d = Double.parseDouble(value);
				recorded_vector_values.add(d + current_offset);
				if (secondary_recorded_values != null) secondary_recorded_values.add(d);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	@Override
	public void removeAbsolutePaths(File f) throws IOException {
		// Do nothing as paths, string delimiters language specific
		// so method needs to be overridden in a parser instance.
	}
	
	@Override
	protected void rootLeafHandler(Object context, Node leaf, PrintWriter fout) {
		if (leaf == null) throw new NullPointerException("Tree leaf is NULL.");

		if (leaf.data != null) {
			boolean inPiecewise = false;
			if (isPiecewise(leaf.data)) inPiecewise = true;
			
			if (!isString_(leaf.data)) leaf.data = getSymbol(leaf.data);
			
			String current_value = "", current_symbol = "_FAKE_FAKE_FAKE_FAKE_";
			if (context instanceof IndividualParameterAssignment) {
				IndividualParameterAssignment ipa = (IndividualParameterAssignment) context;
				String format = "\t%s <- %s\n";
				current_value = String.format(format, z.get(ipa.parameter), leaf.data);
			} else if (isLocalVariable(context)) {
				String format = "\t%s <- %s\n";
				current_value = String.format(format, getSymbol(context), leaf.data);
			} else if (isDerivative(context)) {
				DerivativeVariable dv = (DerivativeVariable) context;
				String format = "\t%s <- %s; %s %s\n", description = dv.getSymbId();
				current_symbol = getSymbol(context);
				current_value = String.format(format, current_symbol, leaf.data, comment_char, description);
			} else {
				String format = " %s ";
				current_value = String.format(format, (String) leaf.data);
			}
			
			if (current_value != null) {
				if(inPiecewise) {
					if (current_symbol != null) current_value = current_value.replaceAll(field_tag, current_symbol) + "\n";
				}
				fout.write(current_value);
			}
		} else
			throw new IllegalStateException("Should be a statement string bound to the root.data element.");
	}
	
	private void setAlgorithm(OptimalDesignStep_ step) {
		if (step == null) return;
		
		String definition = step.getAlgorithm();
		if (definition == null) return;
			
		algorithm = OptimisationAlgorithm.fromValue(definition);
	}
	
	/**
	 * Set the type of FIM matrix generated by PFIM.
	 * @param type FIM Type
	 */
	public void setFIMType(FIMtype type) {
		if (type != null) fim_type = type;
	}
	
	/**
	 * Set the type of FIM matrix generated by PFIM.
	 * @param type FIM Type
	 */
	public void setFIMType(String type) {
		if (type != null) fim_type = FIMtype.valueOf(type.toUpperCase());
	}
	
	/**
	 * Specify the R-generated output filename in the CWD.
	 * @param filename Output Filename
	 */
	public void setOutputFIMFilename(String filename) {
		if (filename != null) outputFIMFilename = filename;
	}
	
	/**
	 * Set the program directory (installation path) for the PFIM R package.
	 * @param dir_path
	 */
	public void setProgramDirectory(String dir_path) {
		if (dir_path != null) {
			programDirectory = dir_path;
			programDirectory = programDirectory.replace("\\", PREFERRED_SEPERATOR);
		}
	}
	
	private void setRecordVectorValues(boolean decision, double offset) {
		record_vector_values = decision;
		current_offset = offset;
	}
	
	/**
	 * Specify the STDOUT filename
	 * @param filename
	 */
	public void setStdoutFilename(String filename) {
		if (filename != null) stdoutFilename = filename;
	}
	
	private void writeAdmissibleSamplingTimes(PrintWriter fout, Observation window, char letter) {
		if (fout == null) return;
		
		List<Double> numbers = ob_sampling_map.get(window);
		if (numbers == null) return;
		if (numbers.isEmpty()) return;
		
		Double min = 0.0;
		Double max = numbers.get(0);
		for (Double value : numbers) {
			if (value == null) continue;
			if (value < min) min = value;
			if (value > max) max = value;
		}
		
		String format = "lower%s<-c(%s)\n"; 
		fout.write(String.format(format, letter, min));
	
		format = "upper%s<-c(%s)\n"; 
		fout.write(String.format(format, letter, max));
	}
	
	private void writeAlgorithmOption(PrintWriter fout) {
		if (fout == null) return;
		String format = "algo.option<-\"%s\"\n";
		fout.write(String.format(format, algorithm));
	}
	
	private void writeAllModelFunctions(File output_dir) throws IOException {
		StructuralBlock sb = lexer.getStrucuturalBlock();	
		String output_filename = getModelFunctionFilename(output_dir.getAbsolutePath(), sb);
		
		PrintWriter mout = new PrintWriter(output_filename);	
		writeModelFunction(mout, sb);	
		mout.close();
	}
	
	private void writeAlpha(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "alpha<-%s\n";
		fout.write(String.format(format, alpha));
	}
	
	private void writeAnalyticalModelFunction(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		String idv = z.get(lexer.getAccessor().getIndependentVariable());
		
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		if (tdb == null) throw new NullPointerException("The trial design block was NULL.");
		
		List<PharmMLRootType> targets = tdb.getDoseTargets();
		int variable_count = 0;
		VariableDefinition dose_variable = null;
		for (PharmMLRootType target : targets) {
			if (isLocalVariable(target)) {
				dose_variable = (VariableDefinition) target;
				variable_count++;
			}
		}
		
		if (variable_count != 1) throw new IllegalStateException("PFIM can only support 1 dose target variable for analytical models (AF).");
		z.addReservedWord(dose_variable.getSymbId(), analytical_model_dose_variable);
		
		String format = "%s <- function(%s,%s,%s) {\n";
		fout.write(String.format(format, "form", idv, param_model_symbol, analytical_model_dose_variable));
		
		writeIndividualParameterAssignments(fout);
		writeLocalVariableAssignments(fout, sb);
		writeAnalyticalReturnStatement(fout, sb);
		
		fout.write("}\n");
	}
	
	private void writeAnalyticalReturnStatement(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		ObservationBlock ob = lexer.getObservationBlocks().get(0);
		if (ob == null) 
		throw new IllegalStateException("No observation model declared in the PharmML so unable to determine variable(s) to export.");
		
		List<SimulationOutput> outputs = ob.getSimulationOutputs();
		if (outputs.isEmpty()) 
		throw new IllegalStateException("The observation model has no declared output variables. Unable to generate model function return statement");
		
		if (outputs.size() == 1) {
			String format = "\treturn(%s)\r\n";
			fout.write(String.format(format, z.get(outputs.get(0))));
		} else {
			throw new UnsupportedOperationException("Multiple output variables in analytical functions not supported yet.");
		}
	} 
	
	/**
	 * Create a windows batch file.
	 * @throws IOException
	 */
	public void writeBatchFile() throws IOException {
		String outFilepath = getBatchFilepath();
		
		PrintWriter fout = new PrintWriter(outFilepath);
		fout.write("@echo off\r\n");
		fout.write("rscript --vanilla call_run.r\r\n");
		fout.close();
	}
	
	private void writeBeta(PrintWriter fout) {
		if (fout == null) return;
		
		Accessor a = lexer.getAccessor();
		ParameterBlock pb = lexer.getParameterBlock();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		if (step == null) throw new NullPointerException("OD step is NULL");
		
		List<String> betas = new ArrayList<String>();
		for (IndividualParameter ip : pb.getIndividualParameters()) {
			String beta = "0.0";
			if (step.isFixed(ip)) {
				Object content = ip.getAssign().getContent();
				if (!isSymbolReference(content)) throw new IllegalStateException("A fixed estimation parameter is not referencing numeric parameter");
				SymbolRef ref = (SymbolRef) content;
				PharmMLRootType element = a.fetchElement(ref);
				if (!isPopulationParameter(element)) throw new IllegalStateException("Model element is not the expected population parameter.");
				PopulationParameter p = (PopulationParameter) element;
				referenced_params_map.put(ip, p);
				ParameterEstimate pe = step.getParameterEstimate(p);
				if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
				beta = parse(ctx, lexer.getStatement(pe.getInitialEstimate())).trim();
			} else {
				StructuredModel model = ip.getStructuredModel();
				if (model == null) throw new NullPointerException("The expected structured model is NULL");
				LinearCovariate lcov = model.getLinearCovariate();
				if (lcov == null) throw new NullPointerException("An expected linear covariate was NULL");
				PopulationValue pop_value = lcov.getPopulationValue();
				Object content = pop_value.getAssign().getContent();
				if (!isSymbolReference(content)) throw new IllegalStateException("A population parameter in a structured model is not referencing numeric parameter");
				SymbolRef ref = (SymbolRef) content;
				PharmMLRootType element = a.fetchElement(ref);
				if (!isPopulationParameter(element)) throw new IllegalStateException("Model element is not the expected population parameter.");
				PopulationParameter p = (PopulationParameter) element;
				referenced_params_map.put(ip, p);
				ParameterEstimate pe = step.getParameterEstimate(p);
				if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
				beta = parse(ctx, lexer.getStatement(pe.getInitialEstimate())).trim();
			}
			betas.add(beta);
		}
		
		String format = "beta<-%s\n";
		fout.write(String.format(format, cat(betas)));
	}
	
	private void writeBetaCovariate(PrintWriter fout) {
		if (fout == null) return;
		String betaCovariateValue = "";
		
		String format = "beta.covariate<-list(%s)\n";
		fout.write(String.format(format, betaCovariateValue));
	}
	
	private void writeBetaCovariateOccassions(PrintWriter fout) {
		if (fout == null) return;
		String betaCovariateOccassionValues = "";
		
		String format = "beta.covariate_occ<-list(%s)\n";
		fout.write(String.format(format, betaCovariateOccassionValues));
	}
	private void writeBetaFixed(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		
		List<String> values = new ArrayList<String>();
		for (IndividualParameter ip : pb.getIndividualParameters()) {
			String value = "T";
			if (referenced_params_map.containsKey(ip)) {
				PopulationParameter p = referenced_params_map.get(ip);
				if (step.getFixedParameter(p) == null) value = "F";
			}
			values.add(value);
		}
		
		String format = "beta.fixed<-%s\n";
		fout.write(String.format(format, cat(values)));
	}
	
	private void writeBounds(PrintWriter fout) {
		if (fout == null) return;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb.isODE()) {
			writeBoundsDefault(fout);
			return;
		}
	}
	
	private void writeBoundsDefault(PrintWriter fout) {
		if (fout == null) return;
		fout.write("boundA<-list(c(0,Inf))\n");
	}
	
	private void writeComputeNNI(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(computeNNI)).trim();
		
		String format = "compute.nni<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeComputeNNIEquivalence(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(computeNNIEquivalence)).trim();
		
		String format = "compute.nni_eq<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeComputePower(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(computePower)).trim();
		
		String format = "compute.power<-%s\n";
		fout.write(String.format(format, value));
	} 
	
	private void writeComputePowerEquivalence(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(computePowerEquivalence)).trim();
		
		String format = "compute.power_eq<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeCondInitExpression(PrintWriter fout) {
		if (fout == null) return;
		
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		StringBuffer s = new StringBuffer("condinit<-expression(");
		int i = 0;
		for (ArmDefinition arm : tdb.getArms()) {
			String dosing_vector = dosing_vector_map.get(arm);
			if (dosing_vector == null) throw new NullPointerException("No dose vector generated from the initial conditions");
			if (i > 0) s.append(",");
			s.append(dosing_vector);
			i++;
		}
		s.append(")\n");
		fout.write(s.toString());
	}
	
	private void writeCondInitIdentical(PrintWriter fout) throws IOException {
		if (fout == null) return;
		
		TreeMaker tm = lexer.getTreeMaker();
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		StructuralBlock sb = lexer.getStrucuturalBlock();
		String default_ic_statement = "0.0";
		
		Set<String> unique_dose_vectors = new HashSet<String>();
		for (ArmDefinition arm : tdb.getArms()) {
			InterventionSequenceRef iref = tdb.getInterventionSequenceRef(arm);
			double dose_start = 0.0;
			if (iref != null) dose_start = iref.start;
			
			dose_start += tdb.getAdministrationStartTime(iref.administration_oid);
			PharmMLRootType target = tdb.getDoseTarget(iref.administration_oid);
			if (target == null) throw new NullPointerException("Administration target is NULL");

			if (sb.isODE()) {
				List<DerivativeVariable> dvs = sb.getStateVariables();
				List<String> ic_statements = new ArrayList<String>();
				for (DerivativeVariable dv : dvs) {
					String ic_statememt = default_ic_statement;
					InitialCondition ic = dv.getInitialCondition();
					if (ic != null) {
						if (ic.getInitialValue() != null) {
							ic_statememt = parse(ctx, tm.newInstance(ic.getInitialValue())).trim();
							if (dose_start == 0.0 && target.equals(dv)) {
								ic_statememt = tdb.getDoseStatement(iref.administration_oid);
							}
						}
					}
					ic_statements.add(ic_statememt);
				}
				
				String dose_vector = cat(ic_statements);
				dosing_vector_map.put(arm, dose_vector);
				unique_dose_vectors.add(dose_vector);
			} else {
				// AF so just use a default empty vector.
				List<String> ic_statements = new ArrayList<String>();
				ic_statements.add("0.0");
				String dose_vector = cat(ic_statements);
				dosing_vector_map.put(arm, dose_vector);
				unique_dose_vectors.add(dose_vector);
			}
		}
		
		String result = "TRUE";
		if (unique_dose_vectors.size() > 1) result = "FALSE";
		else if (unique_dose_vectors.size() == 0) result = "TRUE";
		String format = "condinit.identical<-%s\n";
		fout.write(String.format(format, result));
	}
	
	private void writeCovariateCategory(PrintWriter fout) {
		if (fout == null) return;
		String categoryNameValues = "";
		
		String format = "covariate.category<-list(%s)\n";
		fout.write(String.format(format, categoryNameValues));
	}
	
	private void writeCovariateModel(PrintWriter fout) {
		if (fout == null) return;
		String covariateModelFlag = "F";
		
		String format = "covariate.model<-%s\n";
		fout.write(String.format(format, covariateModelFlag));
	}
	
	private void writeCovariateName(PrintWriter fout) {
		if (fout == null) return;
		String covariateNameValues = "";
		
		String format = "covariate.name<-list(%s)\n";
		fout.write(String.format(format, covariateNameValues));
	}
	
	private void writeCovariateOccassionCategory(PrintWriter fout) {
		if (fout == null) return;
		String covariateOccassionCategoryValues = "";
		
		String format = "covariate_occ.category<-list(%s)\n";
		fout.write(String.format(format, covariateOccassionCategoryValues));
	}
	
	private void writeCovariateOccassionModel(PrintWriter fout) {
		if (fout == null) return;
		String covariateOccassionModel = "F";
		
		String format = "covariate_occ.model<-%s\n";
		fout.write(String.format(format, covariateOccassionModel));
	}
	
	private void writeCovariateOccassionName(PrintWriter fout) {
		if (fout == null) return;
		String covariateOccassionNames = "";
		
		String format = "covariate_occ.name<-list(%s)\n";
		fout.write(String.format(format, covariateOccassionNames));
	}
	
	private void writeCovariateOccassionProportions(PrintWriter fout) {
		if (fout == null) return;
		String covariateOccassionProportions = "";
		
		String format = "covariate_occ.proportions<-list(%s)\n";
		fout.write(String.format(format, covariateOccassionProportions));
	}
	
	private void writeCovariateOccassionSequence(PrintWriter fout) {
		if (fout == null) return;
		String covariateOccassionSequenceValues = "";
		
		String format = "covariate_occ.sequence<-list(%s)\n";
		fout.write(String.format(format, covariateOccassionSequenceValues));
	}
	
	private void writeCovariateProportions(PrintWriter fout) {
		if (fout == null) return;
		String covariateProportionValue = "";
		
		String format = "covariate.proportions<-list(%s)\n";
		fout.write(String.format(format, covariateProportionValue));
	}
	
	private void writeDeltaTime(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "delta.time<-%s\n";
		fout.write(String.format(format, deltaTime));
	}
	
	private void writeDerivatives(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		List<DerivativeVariable> dvs = sb.getStateVariables();
		if (dvs.isEmpty()) return;
		
		for (DerivativeVariable dv : dvs) parse(dv, lexer.getStatement(dv), fout);
		fout.write("\n");
	}
	
	private void writeDose(PrintWriter fout) {
		if (fout == null) return;
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		
		List<ArmDefinition> arms = tdb.getArms();
		if (arms.isEmpty()) return;
		
		int i = 0;
		StringBuffer doses = new StringBuffer();
		for (ArmDefinition arm : arms) {
			if (arm == null) continue;
			if (i > 0) doses.append(",");
			InterventionSequenceRef iref = tdb.getInterventionSequenceRef(arm);
			if (iref == null) throw new NullPointerException("Intervention sequence is NULL");
			String stmt = tdb.getDoseStatement(iref.administration_oid);
			if (stmt == null) throw new NullPointerException("Dose statement is NULL");
			doses.append(stmt);
			i++;
		}
		
		String format = "dose<-c(%s)\n";
		fout.write(String.format(format, doses));
	}
	
	private void writeEquivalenceInterval(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "interval_eq<-c(log(0.8),log(1.25))\n";
		fout.write(format);
	}
	
	private void writeErrorModelDerivedSettings(PrintWriter fout) {
		List<ObservationBlock> obs = lexer.getObservationBlocks();
		int i = 0;
		for (ObservationBlock ob : obs) processErrorModel(i++, ob, fout);
	}
	
	private void writeFedorovWynnOptions(PrintWriter fout) {
		if (fout == null) return;
		
		OptimalDesignLexer od_lexer = (OptimalDesignLexer) lexer;
		OptimalDesignStep_ step = od_lexer.getOptimalDesignStep();
		if (step == null) throw new NullPointerException("Optional Design Step is NULL");
		if (!step.isOptimisation()) return;
		
		TrialDesignBlock2 tdb =  (TrialDesignBlock2) lexer.getTrialDesign();
		List<SingleDesignSpace> spaces = tdb.getDesignSpaces();
		if (spaces.size() < 2) throw new IllegalStateException("Not all design spaces assigned by for trial design optimisation.");
		
		writeNumberOfSamplingWindows(fout, spaces);
		writeSamplingWindows(fout, tdb, spaces);
	}
		
	private void writeFileModel(PrintWriter fout) {
		if (fout == null) return;
		String model_filename = getModelFilename();
		String format = "file.model<-\"%s\"\n";
		fout.write(String.format(format, model_filename));
	}
	
	private void writeFIMType(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "FIM<-\"%s\"\n";
		fout.write(String.format(format, fim_type));
	}
	
	private void writeFixedTimes(PrintWriter fout, char letter) {
		if (fout == null) return;
		String format = "fixed.times%s<-c()\n";
		fout.write(String.format(format, letter));
	}
	
	private void writeGamma(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		
		if (!tdb.hasOccassions()) {
			List<String> values = new ArrayList<String>();
			for (int i = 0; i < pb.getIndividualParameters().size(); i++) values.add("0.0");
				
			String format = "gamma<-diag(%s)\n";
			fout.write(String.format(format, cat(values)));
		}
	}
	
	private void writeGivenPower(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(givenPower)).trim();
		
		String format = "given.Power<-%s\n";
		fout.write(String.format(format, value));
	} 
	
	private void writeGraphInf(PrintWriter fout) {
		if (fout == null) return;
		String format = "graph.infA<-c(0)\n";
		fout.write(format);
	}
	
	private void writeGraphLogical(PrintWriter fout) {
		if (fout == null) return;
		fout.write("graph.logical<-T\n");
	}
	
	private void writeGraphOnly(PrintWriter fout) {
		if (fout == null) return;
		String decision = "FALSE";
		if (lexer.hasPlottingBlock()) decision = "TRUE";
		
		String format = "graph.only<-%s\n";
		fout.write(String.format(format, decision));
	}
	
	private void writeGraphSensiLogical(PrintWriter fout) {
		if (fout == null) return;
		fout.write("graphsensi.logical<-F\n");
	}
	
	private void writeGraphSup(PrintWriter fout) {
		if (fout == null) return;
		
		double max_time = 0.0;
		for (double time : recorded_vector_values) if (time > max_time) max_time = time;
		if (sr != null) {
			if (sr.hasValue(GRAPH_SUPA)) max_time = Double.parseDouble(sr.getValue(GRAPH_SUPA));
		}
		
		String format = "graph.supA<-c(%s)\n";
		fout.write(String.format(format, max_time));
	}
	
	private void writeIdenticalDose(PrintWriter fout) {
		BooleanValue value = null;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb == null) throw new NullPointerException("Structural block is NULL");
		
		if (sb.isODE()) value = new TrueBoolean();
		else {
			TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
			if (tdb.getArmCount() == 1) value = new TrueBoolean();
		}
		
		if (value == null) throw new IllegalStateException("Identical dose flag not specified");
		
		String flag = getSymbol(value);
		
		String format = "dose.identical<-%s\n";
		fout.write(String.format(format, flag));
	}
	
	private void writeIdenticalTimes(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(identicalTimes)).trim();
		
		String format = "identical.times<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeIndividualParameterAssignments(PrintWriter fout) {
		if (fout == null) return;
		ParameterBlockImpl pb = (ParameterBlockImpl) lexer.getParameterBlock();
		for (Object o : pb.getIndividualParameterAssignments()) parse(o, lexer.getStatement(o), fout);
		fout.write("\n");
	}
	
	private void writeLocalVariableAssignments(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		List<VariableDefinition> locals = sb.getLocalVariables();
		if (locals.isEmpty()) return;
		
		for (VariableDefinition local : locals) parse(local, lexer.getStatement(local), fout);
		fout.write("\n");
	}

	private void writeLogLogical(PrintWriter fout) {
		if (fout == null) return;
		String format = "log.logical<-F\n";
		fout.write(format);
	}
	
	private void writeMaximumIterations(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "Max.iter<-%s\n";
		fout.write(String.format(format, maximumIterations));
	}
	
	/**
	 * Generate a call to PFIM.
	 * @param fout Output
	 */
	public void writeModelCall(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "PFIM(model.file=\"stdin.r\")";
		fout.write(format);
	}
	
	private void writeModelForm(PrintWriter fout) {
		if (fout == null) return;
		String form = "AF";
		
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb == null) throw new NullPointerException("Structural Block is NULL"); 
		if (sb.isODE()) form = "DE";
		
		String format = "modelform<-\"%s\"\n";
		fout.write(String.format(format, form));
	}
	
	private void writeModelFunction(PrintWriter fout, StructuralBlock sb) throws IOException {
		if (fout == null) throw new NullPointerException();
		
		writeScriptHeader(fout, lexer.getModelFilename());
		if (sb.isODE()) writeODEModelFunction(fout, sb);
		else writeAnalyticalModelFunction(fout, sb);
	}
	
	
	private void writeNamesDataX(PrintWriter fout) {
		if (fout == null) return;
		fout.write("names.datax<-c(\"Time\")\n");
	}
	
	private void writeNamesDataY(PrintWriter fout) {
		if (fout == null) return;
		String yLabel = "Amount";
		
		Converter c = (Converter) lexer;
		List<VariableDefinition> exported_locals = c.getExportedLocalVariables();
		if (!exported_locals.isEmpty()) yLabel = z.get(exported_locals.get(0));
		
		String format = "names.datay<-c(\"%s\")\n";
		fout.write(String.format(format, yLabel));
	}
	
	private void writeNum(PrintWriter fout) {
		if (fout == null) return;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb.isODE()) {
			fout.write("NUM<-F\n");
		} else {
			fout.write("NUM<-T\n");
		}
	}
	
	private void writeNumberOfOccassions(PrintWriter fout) {
		if (fout == null) return;
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		if (!tdb.hasOccassions()) fout.write("n_occ<-1\n");
	}
	
	private void writeNumberOfResponses(PrintWriter fout) {
		if (fout == null) return;
		Integer nr = lexer.getObservationBlocks().size();
		String format = "nr<-%s\n";
		
		fout.write(String.format(format, nr, comment_char));
	}
	
	private void writeNumberOfSamplingWindows(PrintWriter fout, List<SingleDesignSpace> spaces) {
		if (fout == null || spaces == null) return;
		
		int nWindows = 0;
		for (SingleDesignSpace space : spaces) if(space.getObservationTimes() != null) nWindows++; 
		
		String format = "nwindA<-%s\n";
		fout.write(String.format(format, nWindows));
	}
	
	private void writeODEModelFunction(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		String idv = z.get(lexer.getAccessor().getIndependentVariable());
		
		String format = "%s <- function(%s,%s,%s) {\n";
		fout.write(String.format(format, "formED", idv, state_vector_symbol, param_model_symbol));
		
		writeIndividualParameterAssignments(fout);
		writeLocalVariableAssignments(fout, sb);
		writeDerivatives(fout, sb);
		writeODEModelFunctionReturnStatement(fout, sb);
		
		fout.write("}\n");
	}
	private void writeODEModelFunctionReturnStatement(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		StringBuffer dv_array = new StringBuffer("c(");
		int i = 0;
		for (DerivativeVariable dv : sb.getStateVariables()) {
			if (i > 0) dv_array.append(",");
			dv_array.append(getSymbol(dv));
			i++;
		}
		dv_array.append(")");
		
		Converter c = (Converter) lexer;
		List<VariableDefinition> exported_variables = c.getExportedLocalVariables();
		boolean has_exported_variables = exported_variables.size() > 0;
		StringBuffer exported_variables_v = new StringBuffer("c(");
		i = 0;
		for (VariableDefinition v : exported_variables) {
			if (i > 0) exported_variables_v.append(",");
			exported_variables_v.append(getSymbol(v));
			i++;
		}
		exported_variables_v.append(")");
		
		if (has_exported_variables) {
			String format = "\treturn(list(%s,%s))\n";
			fout.write(String.format(format, dv_array, exported_variables_v));
		} else {
			String format = "\treturn(list(%s))\n";
			fout.write(String.format(format, dv_array));
		}
	}
	
	private void writeOmega(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		if (step == null) throw new NullPointerException("OD step is NULL");
		
		List<String> omegas = new ArrayList<String>();
		for (IndividualParameter ip : pb.getIndividualParameters()) {
			String omega = "0.0";
			
			
			if (!step.isFixed(ip)) {
				ParameterRandomEffect rv = ip.getStructuredModel().getListOfRandomEffects().get(0);
				omega = readOmega(rv);
			}
			
			omegas.add(omega);
		}
		
		String format = "omega<-diag(%s)\n";
		fout.write(String.format(format, cat(omegas)));
	}
	
	
	private void writeOptimisationOptions(PrintWriter fout) {
		if (fout == null) return;
		
		OptimalDesignLexer od_lexer = (OptimalDesignLexer) lexer;
		OptimalDesignStep_ step = od_lexer.getOptimalDesignStep();
		
		if (step == null) throw new NullPointerException("Optional Design Step is NULL");
		if (!step.isOptimisation()) return;
		setAlgorithm(step);
		
		writeIdenticalTimes(fout);
		writeAlgorithmOption(fout);
		
		if (FEDOROV_WYNN.equals(algorithm)) writeFedorovWynnOptions(fout);
		else if (SIMPLEX.equals(algorithm)) writeSimplexOptions(fout);
	}
	
	private void writeOptimiseOnSubjectProportion(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(optProportionsOfsubjects)).trim();
		
		String format = "subjects.opt<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeOption(PrintWriter fout) {
		if (fout == null) return;
		String format = "option<-%s\n";
		fout.write(String.format(format, option));
	}
	
	private void writeOutput(PrintWriter fout) {
		if (fout == null) return;
		String format = "output<-\"%s\"\n";
		String filename = stdoutFilename;
		if (sr != null) if (sr.hasValue(OUTPUT)) filename = sr.getValue(OUTPUT);
		fout.write(String.format(format, filename));
	}
	
	private void writeOutputFIMFilename(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "outputFIM<-\"%s\"\n";
		String filename = outputFIMFilename;
		if (sr != null) if (sr.hasValue(OUTPUT_FIM)) filename = sr.getValue(OUTPUT_FIM);
		fout.write(String.format(format, filename));
	}
	
	private void writeParameterAssociated(PrintWriter fout) {
		if (fout == null) return;
		String parameterAssociated = "";
		
		String format = "parameter.associated<-list(%s)\n";
		fout.write(String.format(format, parameterAssociated));
	}

	private void writeParameterOccassionAssociated(PrintWriter fout) {
		if (fout == null) return;
		String parameterOccassionAssociatedValues = "";
		
		String format = "parameter_occ.associated<-list(%s)\n";
		fout.write(String.format(format, parameterOccassionAssociatedValues));
	}
	
	private void writeParameters(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		if (pb == null) return;
		
		List<IndividualParameter> ips = pb.getIndividualParameters();
		List<String> names = new ArrayList<String>();
		String format = "\"%s\"";
		for (IndividualParameter ip : ips) names.add(String.format(format, z.get(ip)));
		
		format = "parameters<-%s\n";
		fout.write(String.format(format, cat(names)));
	}
	
	private void writePFIMProjectFile() throws IOException {
		if (programDirectory == null) throw new NullPointerException("PFIM Program Directory not specified.");
		String outputFilepath = getPFIMProjectFilepath();
		PrintWriter fout = new PrintWriter(outputFilepath);
		
		writeScriptHeader(fout, null);
		
		// Filter the project template for project specific settings.
		for (int i = 0; i < pfimProjectTemplate.size(); i++) {
			String line = pfimProjectTemplate.get(i);
			if (line == null) continue;
			
			if (line.contains(directory_replacement_label))
				line = line.replace(directory_replacement_label, lexer.getOutputDirectory());
			else if (line.contains(directory_program_replacement_label))
				line = line.replace(directory_program_replacement_label, programDirectory);
			
			fout.write(line + "\n");
		}
		fout.close();
	}
	
	
	@Override
	public void writePreMainBlockElements(PrintWriter fout, File output_dir) throws IOException{
		if (fout == null) return;
		
		writePFIMProjectFile();
		writeScriptHeader(fout, lexer.getModelFilename());
		writeAllModelFunctions(output_dir);
		writeScriptLibraryReferences(fout);
	}
	
	private void writePreviousFIM(PrintWriter fout) {
		if (fout == null) return;
		String format = "previous.FIM<-\"%s\"\n";
		fout.write(String.format(format, previousFIM));
	}
	
	private void writePrintIterations(PrintWriter fout) {
		if (fout == null) return;
		
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(printIterations)).trim();
		
		String format = "iter.print<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeProject(PrintWriter fout) {
		if (fout == null) return;
		String format = "project<-\"%s\"\n";
		String project = lexer.getModelName();
		if (sr != null) if (sr.hasValue(PROJECT)) project = sr.getValue(PROJECT);
		fout.write(String.format(format, project));
	}
	private void writeProt(PrintWriter fout) {
		if (fout == null) return;
		
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		List<Observation> obs = tdb.getObservations();
		if (obs.isEmpty()) return;
		
		Map<ArmDefinition, String> sampling_vectors_raw = new HashMap<ArmDefinition, String>(); 
		for (ArmDefinition arm : tdb.getArms()) {
			if (arm == null) continue;
			Observation ob =  tdb.getObservation(arm);
			if (ob == null) continue;
			if (ob.getObservationTimes() == null) continue;
			
			double start = tdb.getObservationStart(arm);
			setRecordVectorValues(true, start);
			secondary_recorded_values = new ArrayList<Double>();
			String sampling = parse(ctx, lexer.getStatement(ob.getObservationTimes())).trim();
			sampling_vectors_raw.put(arm, sampling);
			ob_sampling_map.put(ob, secondary_recorded_values);
			secondary_recorded_values = null;
		}
		
		List<String> arm_sampling_vectors = new ArrayList<String>();
		for (ArmDefinition arm : tdb.getArms()) {
			if (arm == null) continue;
			if (!sampling_vectors_raw.containsKey(arm)) continue;
			String sampling_vector_raw = sampling_vectors_raw.get(arm);
			
			// Remove the 'R' vector notation.
			sampling_vector_raw = sampling_vector_raw.replace("c(", "");
			sampling_vector_raw = sampling_vector_raw.replace(")", "");
			
			double start = tdb.getObservationStart(arm);
			String [] tokens = sampling_vector_raw.split(",");
			List<Double> values = new ArrayList<Double>(); 
			for (String token : tokens) {
				token = token.trim();
				Double value = Double.parseDouble(token) + start;
				values.add(value);
			}
			
			Collections.sort(values);	
			arm_sampling_vectors.add(cat_(values));
		}
		
		if (arm_sampling_vectors.isEmpty()) return;
		StringBuffer s = new StringBuffer();
		int i = 0;
		for (String arm_sampling_vector : arm_sampling_vectors) {
			if (arm_sampling_vector == null) continue;
			if (i > 0) s.append(",");
			s.append(arm_sampling_vector);
			i++;
		}
		
		String format = "protA<- list(%s)\n";
		fout.write(String.format(format, s));
	} 
	private void writeRelativeConvergenceTolerance(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "Rctol<-%s\n";
		fout.write(String.format(format, relativeConvergenceTolerance));
	}
	private void writeRun(PrintWriter fout) {
		if (fout == null) return;
		OptimalDesignLexer od_lexer = (OptimalDesignLexer) lexer;
		OptimalDesignStep_ step = od_lexer.getOptimalDesignStep();
		
		if (step == null) throw new NullPointerException("Optional Design Step is NULL");
		
		String run = null;
		if (step.isEvaluation()) run = "EVAL";
		else if (step.isOptimisation()) run = "OPT";
		
		if (run == null) throw new IllegalStateException("Unable to determine the PFIM 'RUN' option.");
		
		String format = "run<-\"%s\"\n";
		fout.write(String.format(format, run));
	}
	private void writeSamplingPointsPerSubject(PrintWriter fout, char letter) {
		if (fout == null) return;
		if (recorded_vector_values.isEmpty()) return;
		Double min = recorded_vector_values.get(0);
		Double max = recorded_vector_values.get(0);
		
		for (Double value : recorded_vector_values) {
			if (value < min) min = value;
			if (value > max) max = value;
		}
		
		String format = "nmaxpts%s<-%s\n";
		fout.write(String.format(format, letter, max.intValue()));
		
		format = "nminpts%s<-%s\n";
		fout.write(String.format(format, letter, min.intValue()));
	}
	
	private void writeSamplingWindows(PrintWriter fout, TrialDesignBlock2 tdb, List<SingleDesignSpace> spaces) {
		if (fout == null || tdb == null) return;
		
		List<Observation> windows = new ArrayList<Observation>();
		for (SingleDesignSpace space : spaces) {
			Observation ob = tdb.getObservation(space);
			if (ob == null) continue;
			if (!windows.contains(ob)) windows.add(ob);
		}
		
		if (windows.size() == 0) return;
		for (Observation window : windows) {
			int idx = tdb.getObservationIndex(new OidRef(window));
			if (idx == -1) throw new IllegalStateException("Observation index is -1 (broken window reference in PharmML).");
			char letter = (char) ('A' + idx);
			List<String> sampling_windows = new ArrayList<String>();
			for (SingleDesignSpace space : spaces) {
				if (space.getObservationTimes() != null) {
					Observation ob = tdb.getObservation(space);
					if (window.equals(ob)) {
						String stmt = parse(ctx, lexer.getStatement(space.getObservationTimes())).trim();
						sampling_windows.add(stmt);
					}
				}
			}
			
			StringBuffer sb = new StringBuffer();
			int i = 0;
			for (String sampling_window : sampling_windows) {
				if (i > 0) sb.append(",");
				sb.append(sampling_window);
			}
			
			String format = "sampwin%s<-list(%s)\n";
			fout.write(String.format(format, letter, sb));
			writeFixedTimes(fout, letter);
			
			List<String> sampling_count_limits = new ArrayList<String>();
			recorded_vector_values.clear();
			record_vector_values = true;
			for (SingleDesignSpace space : spaces) {
				if (space.getNumberTimes() != null) {
					Observation ob = tdb.getObservation(space);
					if (window.equals(ob)) {
						String stmt = parse(ctx, lexer.getStatement(space.getNumberTimes())).trim();
						sampling_count_limits.add(stmt);
					}
				}
			}
			record_vector_values = false;
			
			sb = new StringBuffer();
			i = 0;
			for (String sampling_count_limit : sampling_count_limits) {
				if (i > 0) sb.append(",");
				sb.append(sampling_count_limit);
			}
			format = "nsamp%s<-list(%s)\n";
			fout.write(String.format(format, letter, sb));
			writeSamplingPointsPerSubject(fout, letter);
		}
	}
	
	@Override
	protected void writeScriptHeader(PrintWriter fout, String model_file) throws IOException {
		if (fout == null) return;

		String format = "%s Script generated by the PFIM Converter\n";
		fout.write(String.format(format, comment_char));

		format = "%s PFIM 4.0, Copyright (c) Universite Paris Diderot and INSERM (2016)\n";
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
	
	private void writeScriptLibraryReferences(PrintWriter fout) throws IOException {
		if (fout == null) return;
		
		String format = "setwd(\"%s\")\n";
		fout.write(String.format(format, lexer.getOutputDirectory()));
		
		format = "directory<-getwd()\n";
		fout.write(format);
		
		format = "directory.program<-\"%s\"\n";
		fout.write(String.format(format, programDirectory));
		
		format = "source('%s')\n\n";
		fout.write(String.format(format, getPFIMProjectFilepath()));
	}
	
	private void writeSimplexOptions(PrintWriter fout) {
		if (fout == null) return;
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		List<SingleDesignSpace> spaces = tdb.getDesignSpaces();
		if (spaces.isEmpty()) return;
		
		writeOptimiseOnSubjectProportion(fout);
		
		List<Observation> windows = new ArrayList<Observation>();
		for (SingleDesignSpace space : spaces) {
			Observation ob = tdb.getObservation(space);
			if (ob == null) continue;
			if (!windows.contains(ob)) windows.add(ob);
		}

		if (windows.size() == 0) return;
		for (Observation window : windows) {
			int idx = tdb.getObservationIndex(new OidRef(window));
			if (idx == -1) throw new IllegalStateException("Observation index is -1 (broken window reference in PharmML).");
			char letter = (char) ('A' + idx);
			
			if (ob_sampling_map.containsKey(window)) writeAdmissibleSamplingTimes(fout, window, letter);
		}
		
		writeDeltaTime(fout);
		writePrintIterations(fout);
		writeSimplexParameter(fout);
		writeMaximumIterations(fout);
		writeRelativeConvergenceTolerance(fout);
	}
	
	private void writeSimplexParameter(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "simplex.parameter<-%s\n";
		fout.write(String.format(format, simplexParameter));
	}
	
	private void writeSolverSettings(PrintWriter fout) {
		if (fout == null) return;
		String format = "%s<-%s\n";
		fout.write(String.format(format, "RtolEQ", rtol));
		fout.write(String.format(format, "AtolEQ", atol));
		fout.write(String.format(format, "Hmax", "Inf"));
	}
	
	/**
	 * Write a PFIM STDIN file.
	 * @throws IOException
	 */
	public void writeSTDIN() throws IOException {
		String outFilepath = getStdinFilepath();
		
		PrintWriter fout = new PrintWriter(outFilepath);
		writeScriptHeader(fout, lexer.getModelFilename());
		writeProject(fout);
		writeFileModel(fout);
		writeOutput(fout);
		writeOutputFIMFilename(fout);
		writeFIMType(fout);
		writePreviousFIM(fout);
		writeRun(fout);
		writeGraphOnly(fout);
		writeOption(fout);
		writeNumberOfResponses(fout);
		writeModelForm(fout);
		writeIdenticalDose(fout);
		writeDose(fout);
		writeBounds(fout);
		writeNum(fout);
		writeTimeCondInit(fout);
		writeCondInitIdentical(fout);
		writeCondInitExpression(fout);
		writeSolverSettings(fout);
		writeParameters(fout);
		writeBeta(fout);
		writeBetaFixed(fout);
		writeNumberOfOccassions(fout);
		writeTrand(fout);
		writeOmega(fout);
		writeGamma(fout);
		writeErrorModelDerivedSettings(fout);
		writeProt(fout);
		writeSubjects(fout);
		writeSubjectsInput(fout);
		writeCovariateModel(fout);
		writeCovariateName(fout);
		writeCovariateCategory(fout);
		writeCovariateProportions(fout);
		writeParameterAssociated(fout);
		writeBetaCovariate(fout);
		writeCovariateOccassionModel(fout);
		writeCovariateOccassionName(fout);
		writeCovariateOccassionCategory(fout);
		writeCovariateOccassionSequence(fout);
		writeCovariateOccassionProportions(fout);
		writeParameterOccassionAssociated(fout);
		writeBetaCovariateOccassions(fout);
		writeAlpha(fout);
		writeComputePower(fout);
		writeComputeNNI(fout);
		writeEquivalenceInterval(fout);
		writeComputePowerEquivalence(fout);
		writeComputeNNIEquivalence(fout);
		writeGivenPower(fout);
		writeOptimisationOptions(fout);
		writeGraphLogical(fout);
		writeGraphSensiLogical(fout);
		writeNamesDataX(fout);
		writeNamesDataY(fout);
		writeLogLogical(fout);
		writeGraphInf(fout);
		writeGraphSup(fout);
		writeYRange(fout);
		fout.close();
	}

	private void writeSubjects(PrintWriter fout) {
		if (fout == null) return;
		
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		List<ArmDefinition> arms = tdb.getArms();
		if (arms.isEmpty()) return;
		
		List<String> sizes = new ArrayList<String>(); 
		for (ArmDefinition arm : arms) {
			if (arm == null) continue;
			Integer size = tdb.getArmSize(arm.getOid());
			sizes.add(size.toString());
		}
		
		String format = "subjects<-%s\n";
		fout.write(String.format(format, cat(sizes)));
	}
	
	private void writeSubjectsInput(PrintWriter fout) {
		if (fout == null) return;
			
		String format = "subjects.input<-1\n";
		fout.write(format);
	}
	
	private void writeTimeCondInit(PrintWriter fout) {
		if (fout == null) return;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb.isODE()) {
			fout.write("time.condinit<-0\n");
			return;
		}
	}
	
	private void writeTrand(PrintWriter fout) {
		if (fout == null) return;
		String format = "Trand<-%s\n";
		fout.write(String.format(format, trand));
	}
	
	private void writeYRange(PrintWriter fout) {
		if (fout == null) return;
		String format = "y.rangeA<-NULL\n";
		fout.write(format);
	}
}

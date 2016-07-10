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
import static crx.converter.engine.PharmMLTypeChecker.isList;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isLogicalBinaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isLogicalUnaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isPiecewise;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;
import static crx.converter.engine.PharmMLTypeChecker.isReal;
import static crx.converter.engine.PharmMLTypeChecker.isRhs;
import static crx.converter.engine.PharmMLTypeChecker.isSequence;
import static crx.converter.engine.PharmMLTypeChecker.isStandardAssignable;
import static crx.converter.engine.PharmMLTypeChecker.isString;
import static crx.converter.engine.PharmMLTypeChecker.isStructuredError;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isTrue;
import static crx.converter.engine.PharmMLTypeChecker.isUnaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isVariableReference;
import static crx.converter.engine.PharmMLTypeChecker.isVector;
import static crx.converter.engine.scriptlets.BaseScriptlet.function_call_format;
import static crx.converter.engine.scriptlets.BaseScriptlet.tmax;
import static crx.converter.engine.scriptlets.BaseScriptlet.tmin;
import static crx.converter.engine.scriptlets.BaseScriptlet.tspan;
import static eu.ddmore.libpharmml.dom.probonto.ParameterName.STDEV;
import static eu.ddmore.libpharmml.dom.probonto.ParameterName.VAR;
import static inserm.converters.pfim.OptimisationAlgorithm.FEDOROV_WYNN;
import static inserm.converters.pfim.OptimisationAlgorithm.SIMPLEX;
import static inserm.converters.pfim.SettingLabel.GRAPH_LOGICAL;
import static inserm.converters.pfim.SettingLabel.GRAPH_SUPA;
import static inserm.converters.pfim.SettingLabel.OUTPUT;
import static inserm.converters.pfim.SettingLabel.OUTPUT_FIM;
import static inserm.converters.pfim.SettingLabel.PROJECT;
import inserm.converters.pfim.parts.OptimalDesignStepImpl;
import inserm.converters.pfim.parts.ParameterBlockImpl;
import inserm.converters.pfim.parts.TrialDesignBlockImpl;
import inserm.converters.pfim.scriptlets.FactoryImpl;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import javax.xml.bind.JAXBElement;

import org.python.core.PyList;
import org.python.core.PyObject;

import crx.converter.engine.Accessor;
import crx.converter.engine.ConversionDetail_;
import crx.converter.engine.SymbolReader.ModifiedSymbol;
import crx.converter.engine.common.BaseParser;
import crx.converter.engine.common.ElementaryDesign;
import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.engine.common.InterventionSequenceRef;
import crx.converter.engine.common.Protocol;
import crx.converter.engine.common.SimulationOutput;
import crx.converter.engine.scriptlets.Factory;
import crx.converter.engine.scriptlets.Host;
import crx.converter.spi.OptimalDesignLexer;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.blocks.TrialDesignBlock2;
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
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.Sequence;
import eu.ddmore.libpharmml.dom.commontypes.StandardAssignable;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.Symbol;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.SymbolType;
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
import eu.ddmore.libpharmml.dom.trialdesign.SingleDesignSpace;
import eu.ddmore.libpharmml.dom.uncertml.VarRefType;

/**
 * PFIM R-based code generator.
 */
public class Parser extends BaseParser {
	private static final String directory_program_replacement_label = "DIRECTORY_PROGRAM_REPLACE";
	private static final String directory_replacement_label = "DIRECTORY_REPLACE";
	private static final String MODEL_FILESTEM = "model";
	private static final String pfimFIMFilename = "FIM.txt";
	private static final String pfimProjectFilename = "PFIM";
	private static final String pfimStdinFilename = "stdin";
	private static final String pfimStdoutFilename = "stdout.out";
	private static final String PREFERRED_SEPERATOR = "/";
	
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
	private boolean graphLogical = true;
	private boolean graphOnly = false;
	private Host host = null;
	private boolean identicalTimes = true;
	private String leftArrayBracket = null;
	private int maximumIterations = 5000;
	
	/**
	 * Parser Name
	 */
	protected String name = "default";
	
	private OutputFIMFormat option = OutputFIMFormat.BLOCK_DIAGONAL_FIM;
	private boolean optProportionsOfsubjects = false;
	private String output_state_vector_symbol = "yd";
	private String outputFIMFilename = pfimFIMFilename;
	private String param_model_symbol  = null;
	private String pfim_result_variable = "ypkpd";
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
	
	private String cat_(List<String> values) {
		if (values == null) return "";
		StringBuffer s = new StringBuffer();
		int i = 0;
		for (String value : values) {
			if (value == null) continue;
			if (i > 0) s.append(",");
			s.append(value);
			i++;
		}
	
		return s.toString();
	}
	
	private String doArrayAccess(String variableName, Integer idx) {
		String format = "%s%s%s%s";
		
		if (variableName == null || idx == null) throw new NullPointerException("Required array access data is NULL");
		
		return String.format(format, variableName, leftArrayBracket, idx, rightArrayBracket);
	}
	
	/**
	 * Convert a big integer
	 * @param i Big Integer
	 * @return String
	 */
	protected String doBigInteger(BigInteger i) { return i.toString(); }
	
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
	
	protected String doElement(JAXBElement<?> element) {
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
	
	protected String doFalse() { return "FALSE"; }
	
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
	
	/**
	 * Convert a PharmNL integer
	 * @param i Integer
	 * @return String
	 */
	protected String doInt(IntValue i) { return i.getValue().toString(); }
	
	
	
	/**
	 * Convert a Java boolean to R,
	 * @param value Value
	 * @return String
	 */
	protected String doJavaBoolean(Boolean value) {
		if (value) return doTrue();
		else return doFalse();
	}
	
	protected String doLocalVariable(VariableDefinition v) { return z.get(v.getSymbId()); }
	
	private String doLogicalBinaryOperator(LogicBinOp l_b_op) { return getLogicalOperator(l_b_op.getOp()); }
	
	private String doLogicalUnaryOperator(LogicUniOp u_b_op) {return getLogicalOperator(u_b_op.getOp()); }
	
	/**
	 * Generate code for a parameter reference.
	 * @param p Parameter
	 * @return String
	 */
	protected String doParameter(PopulationParameter p) { return z.get(p); }
	
	private String doPiecewise(Piecewise pw) { return unassigned_symbol; }
	
	protected String doPiecewiseForFunctionArgument(FunctionArgument arg) {
		if (!isPiecewiseFunctionArgument(arg)) return unassigned_symbol;
		
		Piecewise pw = arg.getAssign().getPiecewise();
		List<Piece> pieces = pw.getListOfPiece();
		if (pieces.isEmpty()) throw new IllegalStateException("The function argument has no expected piecewise statement.");
		
		TreeMaker tm = lexer.getTreeMaker();
		String format = "(%s)";
		return String.format(format, parse(new Object(), tm.newInstance(pieces.get(0).getCondition())));
	}
	
	protected String doReal(RealValue r) { return Double.toString(r.getValue()); }
	
	protected String doRhs(Rhs eq) {
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
	
	protected String doStringValue(StringValue sv) {
		String format = "'%s'";
		return String.format(format, sv.getValue());
	}
	
	private String doSymbolRef(SymbolRef s) {
		String symbol = s.getSymbIdRef();
		
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
	
	protected String doTrue() { return "TRUE"; }
	
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
			value = value.trim();
			st.append(value);
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
	
	private String getBatchFilepathForUnix() {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + "run.sh";
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
	
	private void initHost() {
		if (host == null) {
			Factory factory = new FactoryImpl();
			host = factory.newInstance();
		}
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
	
	private String parseNumericalListExpression(StandardAssignable expr, double offset) {
		initHost();
		if (expr == null) throw new NullPointerException("PharmML list expression undefined.");
		if (expr.getAssign() == null) throw new NullPointerException("Expression assignment not specified");
		if (!isList(expr)) throw new IllegalStateException("A PharmML expression is not an expected list comprehension");
		
		StringBuffer stmt = new StringBuffer();
		String array_variable = "arr";
		
		stmt.append(String.format("%s = []\n", tspan));
		stmt.append(host.parse(variable(array_variable, null), expr));
		stmt.append(String.format("result = list(set(%s))\n", array_variable));
		stmt.append("result.sort()\n");
		
		host.execute(stmt.toString());
		PyObject values = host.getResult("result");
		Vector sampling_times = vector(parseResultList(values, 0.0));
		
		// Process sampling vector into an R-based format.
		BinaryTree bt = lexer.getTreeMaker().newInstance(sampling_times);
		lexer.addStatement(sampling_times, bt);
		lexer.updateNestedTrees();
		return (parse(ctx, bt).trim());
	}
	
	private List<Double> parseResultList(PyObject result, double offset) {
		if (!(result instanceof PyList)) throw new IllegalStateException("A scriptlet result variable is not the expected list comprehension");
		PyList py_list = (PyList) result;
		List<Double> list = new ArrayList<Double>();
		for(Object o : py_list) {
			Double v = Double.parseDouble(o.toString());
			v += offset;
			list.add(Double.parseDouble(v.toString()));
		}
		
		return list;
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
		} else if (SupportedErrorModel.ADDITIVE_ERROR.equals(model_flag)) {
			for (FunctionArgument arg : call.getListOfFunctionArgument()) {
				if (additive.equals(arg.getSymbId())) {
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
	private String readWidth(ParameterRandomVariable rv) {
		String omega = "0.0";
		if (rv == null) return omega;
		
		// Just assuming probonto and normal usage, nothing else.
		Distribution dist = rv.getDistribution();
		if (dist == null) return omega;
		
		ProbOnto probonto = dist.getProbOnto();
		if (probonto == null) return omega;
		
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStep_ step = c.getOptimalDesignStep();
		
		for (DistributionParameter dp : probonto.getListOfParameter()) {
			if (dp == null) continue;
			ParameterName name = dp.getName();
			if (name == null) continue;
			if (STDEV.equals(name)) {
				Object content = dp.getAssign().getContent();
				if (!isSymbolReference(content)) return omega;
				SymbolRef ref = (SymbolRef) content;
				
				PharmMLRootType element = a.fetchElement(ref);
				if (!isPopulationParameter(element)) return omega;
				PopulationParameter p = (PopulationParameter) element;
				
				ParameterEstimate pe = step.getParameterEstimate(p);
				if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
				omega = parse(ctx, tm.newInstance(pe.getInitialEstimate())).trim();
			} else if (VAR.equals(name)) {
				Object content = dp.getAssign().getContent();
				if (!isSymbolReference(content)) return omega;
				SymbolRef ref = (SymbolRef) content;
				
				PharmMLRootType element = a.fetchElement(ref);
				if (!isPopulationParameter(element)) return omega;
				PopulationParameter p = (PopulationParameter) element;
				
				ParameterEstimate pe = step.getParameterEstimate(p);
				if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
				omega = parse(ctx, tm.newInstance(pe.getInitialEstimate())).trim();
				double variance = Double.parseDouble(parse(ctx, tm.newInstance(pe.getInitialEstimate())).trim());
				Double omega_ = Math.sqrt(variance);
				omega = omega_.toString();
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
	
	/**
	 * Register a settings reader with th PFIM parser.
	 * @param sr_ Settings Reader
	 */
	public void register(SettingReader sr_) {  sr = sr_; }
	
	
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
	
	/**
	 * Specify the STDOUT filename
	 * @param filename
	 */
	public void setStdoutFilename(String filename) {
		if (filename != null) stdoutFilename = filename;
	}
	private String toAlphabetic(int i) {
	    if(i < 0) return "-"+ toAlphabetic(-i-1);
	    
	    int quot = i / 26;
	    int rem = i % 26;
	    char letter = (char)((int) 'A' + rem);
	    if (quot == 0) return ""+letter;
	    else return toAlphabetic(quot-1) + letter;
	}
	
	private VariableDefinition variable(String symbId, Rhs assign) {
		if (symbId == null) return null;
		
		VariableDefinition v = new VariableDefinition();
		v.setSymbId(symbId);
		v.setSymbolType(SymbolType.REAL);
		v.setAssign(assign);
		
		return v;
	}
	
	private Vector vector(List<Double> values) {
		Vector vec = new Vector();
		VectorElements elements = vec.createVectorElements();
		for (Double v : values) elements.createRealValue(v.doubleValue());
		
		return vec;
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
		writeBatchFileWindows();
		writeBatchFileUnix();
	}
	
	private void writeBatchFileUnix() throws IOException {
		String outFilepath = getBatchFilepathForUnix();
		
		PrintWriter fout = new PrintWriter(outFilepath);
		fout.write("rscript --vanilla call_run.r\n");
		fout.close();
	}
	
	private void writeBatchFileWindows() throws IOException {
		String outFilepath = getBatchFilepath();
		 
		PrintWriter fout = new PrintWriter(outFilepath);
		fout.write("@echo off\r\n");
		fout.write("rscript --vanilla call_run.r\r\n");
		fout.close();
	}
	
	private void writeBeta(PrintWriter fout) {
		if (fout == null) return;
		
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
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		List<Protocol> protocols = tdb.getProtocols();
		
		String format = "bound%s<-list(c(0,Inf))\n";
		for (Protocol protocol : protocols)
		fout.write(String.format(format, protocol.getLabel()));
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
		
		writeSamplingWindows(fout, tdb);
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
	
	private void writeFixedTimes(PrintWriter fout, String letter) {
		if (fout == null) return;
		String format = "fixed.times%s<-c()\n";
		fout.write(String.format(format, letter));
	}

	private void writeGammas(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		if (step == null) throw new NullPointerException("OD step is NULL");
		
		List<String> gamma_values = new ArrayList<String>();
		for (IndividualParameter ip : pb.getIndividualParameters()) {
			String gamma = "0.0";
			if (!step.isFixed(ip)) {
				List<ParameterRandomVariable> gammas = step.getGammas(ip);
				if (gammas != null) {
					if (gammas.size() > 1) throw new IllegalStateException("PFIM only allows a single omega to be bound to an Individual parameter");
					if (gammas.isEmpty()) continue;
					gamma = readWidth(gammas.get(0));
				}
			}
			
			gamma_values.add(gamma);
		}
		
		String format = "gamma<-diag(%s)\n";
		fout.write(String.format(format, cat(gamma_values)));
	}
	
	private void writeGivenPower(PrintWriter fout) {
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(givenPower)).trim();
		
		String format = "given.power<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphInf(PrintWriter fout) {
		if (fout == null) return;
		
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		List<Protocol> procs = tdb.getProtocols();
		
		String format = "graph.inf%s<-c(0)\n";
		for (Protocol proc : procs) {
			fout.write(String.format(format, proc.getLabel()));
		}
	}
	
	private void writeGraphLogical(PrintWriter fout) {
		if (sr != null) {
			if (sr.hasValue(GRAPH_LOGICAL)) graphLogical = Boolean.parseBoolean(sr.getValue(GRAPH_LOGICAL)); 
		}
		
		if (fout == null) return;
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(graphLogical)).trim();
		String format = "graph.logical<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphOnly(PrintWriter fout) {
		if (fout == null) return;
		
		TreeMaker tm = lexer.getTreeMaker();
		String value = parse(ctx, tm.newInstance(graphOnly)).trim();
		
		String format = "graph.only<-%s\n";
		fout.write(String.format(format, value));
	}
	
	private void writeGraphSensiLogical(PrintWriter fout) {
		if (fout == null) return;
		fout.write("graphsensi.logical<-F\n");
	}
	
	private void writeGraphSup(PrintWriter fout) {
		if (fout == null) return;
		
		double settings_max_time = 0.0;
		if (sr != null) {
			if (sr.hasValue(GRAPH_SUPA)) settings_max_time = Double.parseDouble(sr.getValue(GRAPH_SUPA));
		}
		
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		List<Protocol> procs = tdb.getProtocols();
		TreeMaker tm = lexer.getTreeMaker();
		
		record_vector_values = true;
		String format = "graph.sup%s<-c(%s)\n";
		for (Protocol proc : procs) {
			if (proc == null) continue;
			recorded_vector_values.clear();
			double max_time = 0.0;
			
			for (ElementaryDesign ed : proc.elementary_designs) {
				if (ed == null) continue;
				for (PharmMLRootType expr : ed.timepoint_expressions) {
					lexer.addStatement(expr, tm.newInstance(expr));
					lexer.updateNestedTrees();
					parse(ctx, lexer.getStatement(expr));
				}
			}
			
			for (double time : recorded_vector_values) if (time > max_time) max_time = time;
			if (settings_max_time > 0.0) max_time = settings_max_time;
			
			fout.write(String.format(format, proc.getLabel(), max_time));
		}
		record_vector_values = false;
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
		
		Converter c = (Converter) lexer;
		OptimalDesignStep_ step = c.getOptimalDesignStep();
		if (step.isOptimisation()) {
			String format = "%s<-PFIM(model.file=\"stdin.%s\")";
			fout.write(String.format(format, pfim_result_variable, script_file_suffix));
		} else {
			String format = "%s<-PFIM(model.file=\"stdin.%s\")";
			fout.write(String.format(format, pfim_result_variable, script_file_suffix));
		}
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
		String default_xLabel = "Time";
		List<String> labels = new ArrayList<String>();
		
		Converter c = (Converter) lexer;
		List<Symbol> outputs = c.getContinuousOutputs();
		String format = "\"%s\"";
		for (Symbol output : outputs) 
			if (output != null) labels.add(String.format(format, default_xLabel));
		
		if (labels.isEmpty()) labels.add(default_xLabel);
		format = "names.datax<-%s\n";
		fout.write(String.format(format, cat(labels)));
	}
	
	private void writeNamesDataY(PrintWriter fout) {
		if (fout == null) return;
		String default_yLabel = "Amount";
		List<String> labels = new ArrayList<String>();
		
		Converter c = (Converter) lexer;
		List<Symbol> outputs = c.getContinuousOutputs();
		String format = "\"%s\"";
		for (Symbol output : outputs) 
			if (output != null) labels.add(String.format(format, output.getSymbId()));
		
		if (labels.isEmpty()) labels.add(default_yLabel);
		format = "names.datay<-%s\n";
		fout.write(String.format(format, cat(labels)));
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
	
	private void writeNWind(PrintWriter fout, Protocol protocol) {
		if (fout == null || protocol == null) return;
		String format = "nwind%s<-%s\n";
		fout.write(String.format(format, protocol.getLabel(), protocol.elementary_designs.size()));
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
	
	private void writeOmegas(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		if (step == null) throw new NullPointerException("OD step is NULL");
		
		List<String> omega_values = new ArrayList<String>();
		for (IndividualParameter ip : pb.getIndividualParameters()) {
			String omega = "0.0";
			
			if (!step.isFixed(ip)) {
				List<ParameterRandomVariable> omegas = step.getOmegas(ip);
				if (omegas != null) {
					if (omegas.size() > 1) throw new IllegalStateException("PFIM only allows a single omega to be bound to an Individual parameter");
					if (omegas.isEmpty()) continue;
					omega = readWidth(omegas.get(0));
				}
			}
			omega_values.add(omega);
		}
		
		String format = "omega<-diag(%s)\n";
		fout.write(String.format(format, cat(omega_values)));
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
	
	@Override
	public void writePreMainBlockElements(PrintWriter fout, File output_dir) throws IOException{
		if (fout == null) return;
		
		writeProjectFile();
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
	
	private void writeProjectFile() throws IOException {
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
	
	private void writeProt(PrintWriter fout) {
		if (fout == null) return;
		
		TreeMaker tm = lexer.getTreeMaker();
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		
		List<Protocol> protocols = tdb.getProtocols();
		if (protocols.isEmpty()) throw new IllegalStateException("Trial design defines no sampling protocols");
		
		String catted_arr = "v";
		String protocolFormat = "prot%s<-list(%s)\n";
		initHost();
		for (Protocol protocol : protocols) {
			if (protocol == null) continue;
			List<String> sampling_stmts = new ArrayList<String>();
			int i = 1;
			for (ElementaryDesign ed : protocol.elementary_designs) {
				if (ed == null) continue;
				
				StringBuffer stmt = new StringBuffer();
				stmt.append(String.format("%s = []\n", catted_arr));
				for (PharmMLRootType expr : ed.timepoint_expressions) {
					if (expr == null) continue;
					if (!isStandardAssignable(expr)) throw new IllegalStateException("Unexpected timepoint data type");
					String array_variable = String.format("arr%s", i++);
					stmt.append(host.parse(variable(array_variable, null), expr));
					stmt.append(String.format("%s = %s + %s\n", catted_arr, catted_arr, array_variable));
				}
				stmt.append(String.format(function_call_format, tmin, MIN.toLowerCase(), catted_arr)); 
				stmt.append(String.format(function_call_format, tmax, MAX.toLowerCase(), catted_arr));
				stmt.append(String.format("result = list(set(%s))\n", catted_arr));
				stmt.append("result.sort()\n");
				
				host.execute(stmt.toString());
				String min_time = host.getResult(tmin).toString();
				String max_time = host.getResult(tmax).toString();
				
				PyObject values = host.getResult("result");
	
				Vector sampling_times = vector(parseResultList(values, ed.start_time_offset));
				BinaryTree bt = tm.newInstance(sampling_times);
				lexer.addStatement(sampling_times, bt);
				lexer.updateNestedTrees();
				sampling_stmts.add(parse(ctx, bt).trim());
				
				double min = Double.parseDouble(min_time);
				double max = Double.parseDouble(max_time);
				
				// For the axes limits for PFIM.
				recorded_vector_values.add(min);
				recorded_vector_values.add(max);
			}
			
			fout.write(String.format(protocolFormat, protocol.getLabel(), cat_(sampling_stmts)));
		}
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
	
	private void writeSamplingPointsLimitsForFW(PrintWriter fout, Protocol protocol, Map<ElementaryDesign, List<SingleDesignSpace>> map) {
		if (fout == null || protocol == null || map == null) return;
		
		recorded_vector_values.clear();
		record_vector_values = true;
		TreeMaker tm = lexer.getTreeMaker();
		
		List<String> sampling_limit_vectors = new ArrayList<String>();
		for (ElementaryDesign ed : protocol.elementary_designs) {
			if (ed == null) continue;
			List<SingleDesignSpace> spaces = map.get(ed);
			if (spaces == null) throw new NullPointerException("The design spaces associated with a protocol are NULL (protocol_oid='" + protocol.block + "')");
			if (spaces.isEmpty()) continue;
			
			SingleDesignSpace sampling_limits = null;
			for (SingleDesignSpace space : spaces) {
				if (space == null) continue;
				if (space.getNumberTimes() != null) {
					sampling_limits = space;
					break;
				}
			}
			if (sampling_limits == null) 
				throw new IllegalStateException("Unable to determine the sampling limits for observation window (oid='" + ed.observation_oid + "')");
			
			BinaryTree bt = tm.newInstance(sampling_limits.getNumberTimes());
			lexer.addStatement(sampling_limits.getNumberTimes(), bt);
			lexer.updateNestedTrees();
			sampling_limit_vectors.add(parse(ctx, bt).trim());
		}
		record_vector_values = false;
		
		String format = "nsamp%s<-list(%s)\n";
		fout.write(String.format(format, protocol.getLabel(), cat_(sampling_limit_vectors)));
		
		if (recorded_vector_values.isEmpty()) 
			throw new IllegalStateException("Sampling points limits for a protocol was empty (protocol_oid='" + protocol.block + "')");
		
		Double min = recorded_vector_values.get(0);
		Double max = recorded_vector_values.get(0);
		for (Double value : recorded_vector_values) {
			if (value < min) min = value;
			if (value > max) max = value;
		}
		
		format = "nmaxpts%s<-%s\n";
		fout.write(String.format(format, protocol.getLabel(), max.intValue()));
		
		format = "nminpts%s<-%s\n";
		fout.write(String.format(format, protocol.getLabel(), min.intValue()));
	}
	
	private void writeSamplingWindows(PrintWriter fout, TrialDesignBlock2 tdb) {
		if (fout == null || tdb == null) return;
		
		List<Protocol> protocols = tdb.getProtocols();
		for (Protocol protocol : protocols) {
			writeNWind(fout, protocol);
			Map<ElementaryDesign, List<SingleDesignSpace>> map = tdb.getDesignSpaces(protocol);
			writeSampWin(fout, protocol, map);
			writeFixedTimes(fout, protocol.getLabel());
			writeSamplingPointsLimitsForFW(fout, protocol, map);
		}
	}
	
	private void writeSampWin(PrintWriter fout, Protocol protocol, Map<ElementaryDesign, List<SingleDesignSpace>> map) {
		if (fout == null || protocol == null || map == null) return;
		
		List<String> sampling_vectors = new ArrayList<String>();
		for (ElementaryDesign ed : protocol.elementary_designs) {
			if (ed == null) continue;
			List<SingleDesignSpace> spaces = map.get(ed);
			if (spaces == null) throw new NullPointerException("The design spaces associated with a protocol are NULL (protocol_oid='" + protocol.block + "')");
			if (spaces.isEmpty()) continue;
			
			String sampling_vector = null;
			for (SingleDesignSpace space : spaces) {
				if (sampling_vector != null) break;	
				if (space == null) continue;
				if (space.getObservationTimes() != null) 
					sampling_vector = parseNumericalListExpression(space.getObservationTimes(), ed.start_time_offset);
			}
			if (sampling_vector == null) throw new NullPointerException("Unable to create a sampling vector from a design space linked to an elementary design");
			sampling_vectors.add(sampling_vector);
		}

		String format = "sampwin%s<-list(%s)\n";
		fout.write(String.format(format, protocol.getLabel(), cat_(sampling_vectors)));
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
		
		writeOptimiseOnSubjectProportion(fout);
		writeSimplexSamplingLimits(fout);
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
	
	private void writeSimplexSamplingLimits(PrintWriter fout) {
		if (fout == null) return;
		
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		List<Protocol> protocols = tdb.getProtocols();
		if (protocols.isEmpty()) throw new IllegalStateException("The model contains no sampling protocols");
		TreeMaker tm = lexer.getTreeMaker();
		
		record_vector_values = true;
		for (Protocol protocol : protocols) {
			List<String> lower_limits = new ArrayList<String>(); 
			List<String> upper_limits = new ArrayList<String>();
			if (protocol == null) continue;
			for (ElementaryDesign ed : protocol.elementary_designs) {
				if (ed == null) continue;
				for (PharmMLRootType expr : ed.timepoint_expressions) {
					if (expr == null) continue;
					if (!isList(expr)) 
						throw new IllegalStateException("A PharmML sampling variable is not an expected list comprehension");
					BinaryTree bt = tm.newInstance(expr);
					lexer.addStatement(expr, bt);
					lexer.updateNestedTrees();
					recorded_vector_values.clear();
					parse(ctx, lexer.getStatement(expr));
					
					Double min = 0.0;
					Double max = 0.0;
					for (Double value : recorded_vector_values) {
						if (value < min) min = value;
						if (value > max) max = value;
					}
					
					min += ed.start_time_offset;
					max += ed.start_time_offset;
					
					lower_limits.add(min.toString());
					upper_limits.add(max.toString());
				}
			}
			
			String format = "lower%s<-c(%s)\n"; 
			fout.write(String.format(format, protocol.getLabel(), cat_(lower_limits)));
			
			format = "upper%s<-c(%s)\n"; 
			fout.write(String.format(format, protocol.getLabel(), cat_(upper_limits)));
		}
		record_vector_values = false;
	}
	
	private void writeSolverSettings(PrintWriter fout) {
		if (fout == null) return;
		String format = "%s<-%s\n";
		fout.write(String.format(format, "RtolEQ", rtol));
		fout.write(String.format(format, "AtolEQ", atol));
		fout.write(String.format(format, "Hmax", "Inf"));
	}
	
	private Accessor a = null;
	private TreeMaker tm = null;
	
	/**
	 * Write a PFIM STDIN file.
	 * @throws IOException
	 */
	public void writeSTDIN() throws IOException {
		String outFilepath = getStdinFilepath();
		
		a = lexer.getAccessor();
		tm = lexer.getTreeMaker();
		
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
		writeOmegas(fout);
		writeGammas(fout);
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
		
		TrialDesignBlock2 tdb = (TrialDesignBlock2) lexer.getTrialDesign();
		String format = "y.range%s<-NULL\n";
		for (Protocol proc : tdb.getProtocols()) {
			fout.write(String.format(format, proc.getLabel()));
		}
	}
}

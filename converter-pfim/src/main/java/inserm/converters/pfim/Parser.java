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
import inserm.converters.pfim.parts.OptimalDesignStepImpl;
import inserm.converters.pfim.parts.ParameterBlockImpl;
import inserm.converters.pfim.parts.TrialDesignBlockImpl;
import inserm.converters.pfim.parts.TrialDesignBlockImpl.InterventionSequenceRef;

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
import crx.converter.engine.FixedParameter;
import crx.converter.engine.SymbolReader.ModifiedSymbol;
import crx.converter.engine.common.BaseParser;
import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.spi.OptimalDesignLexer;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.ParameterBlock;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.blocks.TrialDesignBlock;
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
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError.ResidualError;
import eu.ddmore.libpharmml.dom.modellingsteps.FIMtype;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;
import eu.ddmore.libpharmml.dom.probonto.DistributionParameter;
import eu.ddmore.libpharmml.dom.probonto.ParameterName;
import eu.ddmore.libpharmml.dom.probonto.ProbOnto;
import eu.ddmore.libpharmml.dom.trialdesign.ArmDefinition;
import eu.ddmore.libpharmml.dom.trialdesign.Observation;
import eu.ddmore.libpharmml.dom.uncertml.VarRefType;

/**
 * PFIM R-based code generator.
 */
public class Parser extends BaseParser {
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
	
	private Double atol = 1E-08;
	private Object ctx = new Object();
	private double current_offset = 0.0;
	private Map<ArmDefinition, String> dosing_vector_map = new HashMap<ArmDefinition, String>();
	private FIMtype fim_type = FIMtype.P;
	private Double hmax = 1E-08;
	private String leftArrayBracket = null;
	private OutputFIMFormat option = OutputFIMFormat.BLOCK_DIAGONAL_FIM;
	private String output_state_vector_symbol = "yd";
	private String outputFIMFilename = pfimFIMFilename;
	private String param_model_symbol  = null;
	private List<String> pfimProjectTemplate = new ArrayList<String>();
	private String previousFIM = "";
	private String programDirectory = ".";
	private Properties props = null;
	private boolean record_vector_values = false;
	private List<Double> recorded_vector_values = new ArrayList<Double>();
	private Map<IndividualParameter, PopulationParameter> referenced_params_map = new HashMap<IndividualParameter, PopulationParameter>();
	private String rightArrayBracket = null;
	private Double rtol = 1E-08;
	private String state_vector_symbol = null;
	private String stdoutFilename = pfimStdoutFilename;
	private RandomEffectModelOption trand = RandomEffectModelOption.EXPONENTIAL;
	private boolean writtenSTDIN = false;
	
	private boolean wrotePFIM_R = false;
	
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
		script_file_suffix = "r";
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
	
	protected String doFixedParameter(FixedParameter fp) { return unassigned_symbol; }
	
	protected String doFunctionCall(FunctionCallType call) {
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
	
	protected String doIndependentVariable(IndependentVariable v) {
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
	
	@Override
	public String getModelFunctionFilename(String output_dir, StructuralBlock sb) {
		return String.format("%s%s%s.%s", output_dir, File.separator, MODEL_FILESTEM, script_file_suffix);
	}
	
	private String getPFIMProjectFilepath() throws IOException {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + pfimProjectFilename + "." + script_file_suffix;
	}
	
	private String getStdinFilepath() {
		String cwd = lexer.getOutputDirectory();
		return cwd + PREFERRED_SEPERATOR + pfimStdinFilename + "." + script_file_suffix;
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
	
	private boolean is_parameter_scope(ParameterRandomVariable eps) {
		if (eps == null) return false;
		
		LevelReference lref = eps.getListOfVariabilityReference().get(0);
		if (lref != null) {
			VariabilityBlock vb = lexer.getVariabilityBlock(lref.getSymbRef());
			if (vb != null) return vb.isParameterVariability();
		}
		
		return false;
	}
	
	private boolean is_residual_scope(ParameterRandomVariable eps) {
		if (eps == null) return false;
		
		LevelReference lref = eps.getListOfVariabilityReference().get(0);
		if (lref != null) {
			VariabilityBlock vb = lexer.getVariabilityBlock(lref.getSymbRef());
			if (vb != null) return vb.isResidualError();
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

		Accessor a = lexer.getAccessor();
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
		if (!SupportedErrorModel.contains(func_name)) 
			throw new UnsupportedOperationException("Error model not supported (type='" + func_name + "' )");

		String inter = "0.0", slope = "0.0";
		
		OptimalDesignLexer c = (OptimalDesignLexer) lexer;
		OptimalDesignStepImpl step = (OptimalDesignStepImpl) c.getOptimalDesignStep();
		SupportedErrorModel model_flag = SupportedErrorModel.fromValue(func_name);
		if (SupportedErrorModel.PROPORTIONAL.equals(model_flag)) {
			inter = "0.0";
			ResidualError re = serr.getResidualError();
			if (re != null) {
				PharmMLRootType element = a.fetchElement(re.getSymbRef());
				if (!isRandomVariable(element)) throw new IllegalStateException("A model element referenced in a residual error is not a random variable");
				ParameterRandomVariable eps = (ParameterRandomVariable) element;
				if (!is_residual_scope(eps)) throw new IllegalStateException("Random variable does not have residual scope (name='" + eps.getSymbId() +  "')");

				Distribution dist = eps.getDistribution();
				if (dist != null) {
					ProbOnto probo = dist.getProbOnto();
					if (probo != null) {
						for (DistributionParameter dp : probo.getListOfParameter()) {
							if (dp == null) continue;
							ParameterName name = dp.getName();
							if (name == null) continue;
							if (ParameterName.STDEV.equals(name)) {
								content = dp.getAssign().getContent();
								if (!isSymbolReference(content)) continue;
								SymbolRef ref = (SymbolRef) content;
								element = a.fetchElement(ref);
								if (!isPopulationParameter(element)) continue;
								PopulationParameter p = (PopulationParameter) element;

								ParameterEstimate pe = step.getParameterEstimate(p);
								if (pe == null) throw new NullPointerException("The expected parameter estimate is NULL");
								slope = parse(ctx, lexer.getStatement(pe.getInitialEstimate())).trim();
							}
						}
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
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
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
	
	/**
	 * Set the type of FIM matrix generated by PFIM.
	 * @param type FIM Type
	 */
	public void setFIMType(String type) {
		if (type == null) return;
		fim_type = FIMtype.valueOf(type.toUpperCase());
	}
	
	/**
	 * Specify the R-generated output filename in the CWD.
	 * @param filename Output Filename
	 */
	public void setOutputFIMFilename(String filename) {
		if (filename != null) outputFIMFilename = filename;
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
	
	private void writeAllModelFunctions(File output_dir) throws IOException {
		StructuralBlock sb = lexer.getStrucuturalBlock();	
		String output_filename = getModelFunctionFilename(output_dir.getAbsolutePath(), sb);
		
		PrintWriter mout = new PrintWriter(output_filename);	
		writeModelFunction(mout, sb);	
		mout.close();
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

	private void writeCondInitExpression(PrintWriter fout) {
		if (fout == null) return;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		
		if (sb.isODE()) {
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
		} else 
			throw new UnsupportedOperationException("AF not supported yet.");
	}
	
	private void writeCondInitIdentical(PrintWriter fout) throws IOException {
		if (fout == null) return;
		
		TreeMaker tm = lexer.getTreeMaker();
		TrialDesignBlockImpl tdb = (TrialDesignBlockImpl) lexer.getTrialDesign();
		StructuralBlock sb = (StructuralBlock) lexer.getStrucuturalBlock();
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
			} else
				throw new UnsupportedOperationException("AF not supported yet.");
		}
		
		String result = "TRUE";
		if (unique_dose_vectors.size() > 1) result = "FALSE";
		String format = "condinit.identical<-%s\n";
		fout.write(String.format(format, result));
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
	
	private void writeErrorModelDerivedSettings(PrintWriter fout) {
		List<ObservationBlock> obs = lexer.getObservationBlocks();
		int i = 0;
		for (ObservationBlock ob : obs) processErrorModel(i++, ob, fout);
	}
	
	private void writeFIMType(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "FIM<-\"%s\"\n";
		fout.write(String.format(format, fim_type));
	}
	
	private void writeGamma(PrintWriter fout) {
		if (fout == null) return;
		
		ParameterBlock pb = lexer.getParameterBlock();
		TrialDesignBlock tdb = lexer.getTrialDesign();
		
		if (!tdb.hasOccassions()) {
			List<String> values = new ArrayList<String>();
			for (int i = 0; i < pb.getIndividualParameters().size(); i++) values.add("0.0");
				
			String format = "gamma<-diag(%s)\n";
			fout.write(String.format(format, cat(values)));
		}
	}
	
	private void writeGraphOnly(PrintWriter fout) {
		if (fout == null) return;
		String decision = "FALSE";
		if (lexer.hasPlottingBlock()) decision = "TRUE";
		
		String format = "graph.only<-%s\n";
		fout.write(String.format(format, decision));
	}
	
	private void writeIdenticalDose(PrintWriter fout) {
		BooleanValue value = null;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb == null) throw new NullPointerException("Structural block is NULL");
		
		if (sb.isODE()) value = new TrueBoolean();
		
		if (value == null) throw new IllegalStateException("Identical dose flag not specified");
		
		String flag = getSymbol(value);
		
		String format = "dose.identical<-%s\n";
		fout.write(String.format(format, flag));
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
	
	private void writeModelFile(PrintWriter fout) {
		if (fout == null) return;
		String model_filename = MODEL_FILESTEM + "." + script_file_suffix;
		String format = "file.model<-\"%s\"\n";
		fout.write(String.format(format, model_filename));
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
	}
	
	private void writeNum(PrintWriter fout) {
		if (fout == null) return;
		StructuralBlock sb = lexer.getStrucuturalBlock();
		if (sb.isODE()) {
			fout.write("NUM<-F\n");
			return;
		}
	} 
	
	private void writeNumberOfOccassions(PrintWriter fout) {
		if (fout == null) return;
		TrialDesignBlock tdb = lexer.getTrialDesign();
		if (!tdb.hasOccassions()) fout.write("n_occ<-1\n");
	}
	
	private void writeNumberOfResponses(PrintWriter fout) {
		if (fout == null) return;
		Integer nr = lexer.getObservationBlocks().size();
		String format = "nr<-%s\n";
		
		fout.write(String.format(format, nr, comment_char));
	}
	
	private void writeODEModelFunction(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		String idv = z.get(lexer.getAccessor().getIndependentVariable()).toLowerCase();
		
		String format = "%s <- function(%s,%s,%s) {\n";
		fout.write(String.format(format, MODEL_FILESTEM, idv, state_vector_symbol, param_model_symbol));
		
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
	
	private void writeOption(PrintWriter fout) {
		if (fout == null) return;
		String format = "option<-%s\n";
		fout.write(String.format(format, option));
	}
	
	private void writeOutputFilename(PrintWriter fout) {
		if (fout == null) return;
		String format = "output<-\"%s\"\n";
		fout.write(String.format(format, stdoutFilename));
	}
	
	private void writeOutputFIMFilename(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "outputFIM<-\"%s\"\n";
		fout.write(String.format(format, outputFIMFilename));
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
		if (wrotePFIM_R) return;
		
		String outputFilepath = getPFIMProjectFilepath();
		PrintWriter fout = new PrintWriter(outputFilepath);
		
		writeScriptHeader(fout, null);
		
		// Filter the project template for project specific settings.
		for (int i = 0; i < pfimProjectTemplate.size(); i++) {
			String line = pfimProjectTemplate.get(i);
			if (line == null) continue;
			fout.write(line + "\n");
		}
		fout.close();
		
		wrotePFIM_R = true;
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
	
	private void writeProjectName(PrintWriter fout) {
		if (fout == null) return;
		String format = "project<-\"%s\"\n";
		fout.write(String.format(format, lexer.getModelName()));
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
			String sampling = parse(ctx, lexer.getStatement(ob.getObservationTimes())).trim();
			sampling_vectors_raw.put(arm, sampling);
			
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
		String format = "\nsource('%s')\n";
		fout.write(String.format(format, getPFIMProjectFilepath()));
		fout.write("PFIM()\n");
	}
	
	private void writeSolverSettings(PrintWriter fout) {
		if (fout == null) return;
		String format = "%s<-%s\n";
		fout.write(String.format(format, "RtolEQ", rtol));
		fout.write(String.format(format, "AtolEQ", atol));
		fout.write(String.format(format, "Hmax", hmax));
	}
	
	/**
	 * Write a PFIM STDIN file.
	 * @throws IOException
	 */
	public void writeSTDIN() throws IOException {
		if (writtenSTDIN) return;
		
		String outFilepath = getStdinFilepath();
		
		PrintWriter fout = new PrintWriter(outFilepath);
		writeScriptHeader(fout, lexer.getModelFilename());
		writeProjectName(fout);
		writeModelFile(fout);
		writeOutputFilename(fout);
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
		
		fout.close();
		
		writtenSTDIN = true;
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
}

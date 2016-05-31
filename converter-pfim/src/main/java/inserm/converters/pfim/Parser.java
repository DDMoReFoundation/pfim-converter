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
import static crx.converter.engine.PharmMLTypeChecker.isReal;
import static crx.converter.engine.PharmMLTypeChecker.isRhs;
import static crx.converter.engine.PharmMLTypeChecker.isSequence;
import static crx.converter.engine.PharmMLTypeChecker.isString;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isTrue;
import static crx.converter.engine.PharmMLTypeChecker.isUnaryOperation;
import static crx.converter.engine.PharmMLTypeChecker.isVariableReference;
import static crx.converter.engine.PharmMLTypeChecker.isVector;
import inserm.converters.pfim.parts.ParameterBlockImpl;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Properties;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.Accessor;
import crx.converter.engine.ConversionDetail_;
import crx.converter.engine.FixedParameter;
import crx.converter.engine.SymbolReader.ModifiedSymbol;
import crx.converter.engine.common.BaseParser;
import crx.converter.engine.common.IndividualParameterAssignment;
import crx.converter.spi.OptimalDesignLexer;
import crx.converter.spi.blocks.StructuralBlock;
import crx.converter.spi.steps.OptimalDesignStep_;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.Node;
import crx.converter.tree.TreeMaker;
import eu.ddmore.convertertoolbox.api.response.ConversionDetail;
import eu.ddmore.libpharmml.dom.IndependentVariable;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.IntValue;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.Sequence;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
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
import eu.ddmore.libpharmml.dom.modeldefn.PopulationParameter;
import eu.ddmore.libpharmml.dom.modellingsteps.FIMtype;
import eu.ddmore.libpharmml.dom.modellingsteps.ParameterEstimate;
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
	
	private FIMtype fim_type = FIMtype.P;
	private String leftArrayBracket = null;
	private OutputFIMFormat option = OutputFIMFormat.BLOCK_DIAGONAL_FIM;
	private String output_state_vector_symbol = "yd";
	private String outputFIMFilename = pfimFIMFilename;
	private String param_model_symbol  = null;
	private List<String> pfimProjectTemplate = new ArrayList<String>();
	private String previousFIM = "";
	private String programDirectory = ".";
	private Properties props = null;
	private String rightArrayBracket = null;
	private String state_vector_symbol = null;
	private String stdoutFilename = pfimStdoutFilename;
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
		
		ArrayList<String> values = new ArrayList<String>();
		
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
	
	private void writeDerivatives(PrintWriter fout, StructuralBlock sb) {
		if (fout == null || sb == null) return;
		
		List<DerivativeVariable> dvs = sb.getStateVariables();
		if (dvs.isEmpty()) return;
		
		for (DerivativeVariable dv : dvs) parse(dv, lexer.getStatement(dv), fout);
		fout.write("\n");
	}
	
	private void writeFIMType(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "FIM<-\"%s\"\n";
		fout.write(String.format(format, fim_type));
	}
	
	private void writeGraphOnly(PrintWriter fout) {
		String decision = "FALSE";
		if (lexer.hasPlottingBlock()) decision = "TRUE";
		
		String format = "graph.only<-%s\n";
		fout.write(String.format(format, decision));
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
		String model_filename = MODEL_FILESTEM + "." + script_file_suffix;
		String format = "file.model<-\"%s\" %s Model Filename\n";
		fout.write(String.format(format, model_filename, comment_char));
	}
	
	private void writeModelForm(PrintWriter fout) {
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
	
	private void writeNumberOfResponses(PrintWriter fout) {
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
	
	private void writeOption(PrintWriter fout) {
		String format = "option<-%s\n";
		fout.write(String.format(format, option));
	}
	
	private void writeOutputFilename(PrintWriter fout) {
		String format = "output<-\"%s\" %s Output Filename\n";
		fout.write(String.format(format, stdoutFilename, comment_char));
	}

	private void writeOutputFIMFilename(PrintWriter fout) {
		if (fout == null) return;
		
		String format = "outputFIM<-\"%s\"\n";
		fout.write(String.format(format, outputFIMFilename));
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
		String format = "previous.FIM<-\"%s\"\n";
		fout.write(String.format(format, previousFIM));
	}
	
	private void writeProjectName(PrintWriter fout) {
		String format = "project<-\"%s\" %s Project Name\n";
		fout.write(String.format(format, lexer.getModelName(), comment_char));
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
		String format = "\nsource('%s')\n";
		fout.write(String.format(format, this.getPFIMProjectFilepath()));
		fout.write("PFIM()\n");
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
		
		fout.close();
		
		writtenSTDIN = true;
	}
}

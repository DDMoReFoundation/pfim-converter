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

package inserm.converters.pfim.scriptlets;

import static crx.converter.engine.PharmMLTypeChecker.isDerivative;
import static crx.converter.engine.PharmMLTypeChecker.isFalse;
import static crx.converter.engine.PharmMLTypeChecker.isFunction;
import static crx.converter.engine.PharmMLTypeChecker.isInt;
import static crx.converter.engine.PharmMLTypeChecker.isJAXBElement;
import static crx.converter.engine.PharmMLTypeChecker.isLocalVariable;
import static crx.converter.engine.PharmMLTypeChecker.isPiece;
import static crx.converter.engine.PharmMLTypeChecker.isPiecewise;
import static crx.converter.engine.PharmMLTypeChecker.isPopulationParameter;
import static crx.converter.engine.PharmMLTypeChecker.isReal;
import static crx.converter.engine.PharmMLTypeChecker.isRhs;
import static crx.converter.engine.PharmMLTypeChecker.isSequence;
import static crx.converter.engine.PharmMLTypeChecker.isString;
import static crx.converter.engine.PharmMLTypeChecker.isSymbolReference;
import static crx.converter.engine.PharmMLTypeChecker.isTrue;
import static crx.converter.engine.PharmMLTypeChecker.isVector;
import static crx.converter.engine.scriptlets.PyFunctionName.LINSPACE;

import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.ConversionDetail_;
import crx.converter.engine.scriptlets.PyFunctionName;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.Node;
import eu.ddmore.convertertoolbox.api.response.ConversionDetail;
import eu.ddmore.libpharmml.dom.commontypes.IntValue;
import eu.ddmore.libpharmml.dom.commontypes.RealValue;
import eu.ddmore.libpharmml.dom.commontypes.Rhs;
import eu.ddmore.libpharmml.dom.commontypes.Sequence;
import eu.ddmore.libpharmml.dom.commontypes.StringValue;
import eu.ddmore.libpharmml.dom.commontypes.SymbolRef;
import eu.ddmore.libpharmml.dom.commontypes.VariableDefinition;
import eu.ddmore.libpharmml.dom.commontypes.Vector;
import eu.ddmore.libpharmml.dom.commontypes.VectorElements;
import eu.ddmore.libpharmml.dom.commontypes.VectorValue;

/**
 * Mini-python parser to generate statements from PharmML inputs.<br/>
 * List processing in Python is a lot simpler than messing about with python.
 */
public class ScriptLetParser extends inserm.converters.pfim.Parser {
	/**
	 * Constructor
	 */
	public ScriptLetParser() throws IOException {
		super();
		
		comment_char = "#";
		name = "ScriptletHost";
		setReferenceClass(getClass());
		init();
		useConversativeBinaryOperatorFormats();
	}
	
	private String createLinspace() {
		StringBuffer stmt = new StringBuffer();
		
		stmt.append("def linspace(a, b, n=100):\n");
		stmt.append(indent(1) + "if n < 2:\n");
		stmt.append(indent(2) + "return(b)\n");
		stmt.append(indent(1) + "diff = (float(b) - a)/(n - 1)\n");
		stmt.append(indent(1) + "return ([diff * i + a  for i in range(n)])\n\n");
		    
		return stmt.toString();
	}
	
	String createPyFunction(PyFunctionName func) {
		if (LINSPACE.equals(func)) return createLinspace();
		else throw new UnsupportedOperationException("Function Name='" + func + "'");
	}
	
	@Override
	protected String doLocalVariable(VariableDefinition v) { return v.getSymbId(); }
	 
	private String  doSequence(Sequence seq) {
		String symbol = unassigned_symbol;
		
		Rhs begin = seq.getBegin();
		Rhs end = seq.getEnd();
		Rhs repetitions = seq.getRepetitions();
		Rhs stepSize = seq.getStepSize();
		
		BinaryTree btBegin = null, btEnd = null, btRepetitions = null, btStepSize = null;
		if (begin == null) throw new IllegalStateException("The required Sequence.begin field is not assigned.");
		
		if (begin != null) if (lexer.hasStatement(begin)) btBegin = lexer.getStatement(begin);
		if (end != null) if (lexer.hasStatement(end)) btEnd = lexer.getStatement(end);
		if (repetitions != null) if (lexer.hasStatement(repetitions)) btRepetitions = lexer.getStatement(repetitions);
		if (stepSize != null) if (lexer.hasStatement(stepSize)) btStepSize = lexer.getStatement(stepSize);
		
		String strBegin = null, strEnd = null, strRepetitions = null, strStepSize = null;
		
		if (btBegin != null) strBegin = stripOuterBrackets(parse(seq, btBegin)).trim();
		if (btEnd != null) strEnd = stripOuterBrackets(parse(seq, btEnd)).trim();
		if (btRepetitions != null) strRepetitions = stripOuterBrackets(parse(seq, btRepetitions)).trim();
		if (btStepSize != null) strStepSize = stripOuterBrackets(parse(seq, btStepSize)).trim();
		
		// Default value in case the conditional logic fails or the PharmML spec changes.
		symbol = "[0.0, 0.0]";
		if (strBegin != null && strEnd != null && strStepSize != null) {
			String format = "range(%s,%s,%s)";
			symbol = String.format(format, strBegin, strEnd, strStepSize);
		} else if (strBegin != null && strEnd != null && strRepetitions != null) {
			String format = "linspace(%s,%s,%s)";
			symbol = String.format(format, strBegin, strEnd, strRepetitions);
		} else if (strBegin != null && strStepSize != null && strRepetitions != null) {
			String format = "[%s:%s:(%s*%s)]";
			symbol = String.format(format, strBegin, strStepSize, strStepSize, strRepetitions);
		} else if (strBegin != null && strEnd != null) {
			String format = "[%s %s]";
			symbol = String.format(format, strBegin, strEnd);
		}
		
		return symbol;
	}
	
	private String doSymbolRef(SymbolRef s) { return s.getSymbIdRef(); }
	
	private String doVector(Vector v) {
		String symbol = unassigned_symbol;
		List<String> values = new ArrayList<String>();
		
		VectorElements elements = v.getVectorElements();
		if (elements == null) throw new NullPointerException("Vector elements are NULL.");
		
		for (VectorValue value : elements.getListOfElements()) {
			if (value == null) continue;
			if (lexer.hasStatement(value)) values.add(stripOuterBrackets(parse(v, lexer.getStatement(value))));
		}
		
		StringBuilder stmt = new StringBuilder();
		stmt.append("[");
		int i = 0;
		for (String element : values) {
			if (i > 0) stmt.append(",");
			stmt.append(element);
			i++;
		}
		stmt.append("]");
		symbol = stmt.toString();
		
		return symbol;
	}
	
	@Override
	public String getSymbol(Object o) {
		String symbol = unassigned_symbol;
		
		if (isSymbolReference(o)) symbol = doSymbolRef((SymbolRef) o);
		else if (isLocalVariable(o)) symbol = doLocalVariable((VariableDefinition) o);
		
		else if (isReal(o)) symbol = doReal((RealValue) o);
		else if (isFalse(o)) symbol = doFalse();
		else if (isTrue(o)) symbol = doTrue();
		else if (isString(o)) symbol = doStringValue((StringValue) o);
		else if (isInt(o)) symbol = doInt((IntValue) o);	
		else if (isSequence(o)) symbol = doSequence((Sequence) o); 
		else if (isVector(o)) symbol = doVector((Vector) o);
	    else if (isJAXBElement(o)) symbol = doElement((JAXBElement<?>) o);	
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
	
	private String indent(int i) {
		String str = "";
		for (int level = 0; level < i; level++) str += "   ";
		return str;
	}
		
	@Override
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
		lexer.setIndexFromZero(true);
	}

	@Override
	protected void rootLeafHandler(Object context, Node leaf, PrintWriter fout) {
		if (leaf == null) throw new NullPointerException("Tree leaf is NULL.");
		
		if (leaf.data != null) {
			boolean inPiecewise = false;
			if (isPiecewise(leaf.data)) inPiecewise = true;
			
			if (!isString_(leaf.data)) leaf.data = getSymbol(leaf.data);
			String current_value = "", current_symbol = "_FAKE_FAKE_FAKE_FAKE_";
			if (isDerivative(context) || isPopulationParameter(context) || isLocalVariable(context)) {
				String format = "%s = %s\n";
				current_symbol = getSymbol(context);
				current_value = String.format(format, current_symbol, leaf.data);
			} else if (isFunction(context) || isSequence(context) || isVector(context)) {
				String format = "(%s) ";
				current_value = String.format(format, (String) leaf.data);
			} else if (isPiece(context)) { 
				String format = "%s ";
				current_value = String.format(format, (String) leaf.data);
			}  else {
				String format = " %s ";
				current_value = String.format(format, (String) leaf.data);
			} 
			
			if (current_value != null) {
				if (inPiecewise) {
					if (current_symbol != null) current_value = current_value.replaceAll(field_tag, current_symbol) + "\n";
				}
				fout.write(current_value);
			}
		} else
			throw new IllegalStateException("Should be a statement string bound to the root.data element.");
	}
	
	@Override
	public String toString() {
		String format = "%s=%s";
		StringBuffer sb = new StringBuffer();
		sb.append(String.format(format, "name", name));
		return sb.toString();
	}
}

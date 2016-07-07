/*******************************************************************************
 * Copyright (C) 2016 Cyprotex Discovery Ltd - All rights reserved.
 ******************************************************************************/

package inserm.converters.pfim.scriptlets;

import inserm.converters.pfim.Converter;

import java.io.IOException;

import org.python.core.PyObject;
import org.python.core.PySystemState;
import org.python.util.PythonInterpreter;

import crx.converter.engine.scriptlets.Host;
import crx.converter.engine.scriptlets.PyFunctionName;
import crx.converter.spi.IParser;
import crx.converter.tree.BaseTreeMaker;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.NestedTreeRef;
import crx.converter.tree.TreeMaker;
import eu.ddmore.convertertoolbox.domain.LanguageVersionImpl;
import eu.ddmore.convertertoolbox.domain.VersionImpl;

/**
 * A mostly empty class that wraps calls to the Jython scripting run-time.<br/>
 * This class implements ILexer only to provide enough methods so that
 * it can link to the expression engine of BaseParser.<br/>
 * A scriptlet is a small piece of code to run a specific set of operations on data derived
 * from a PharmML model.<br/>
 * As such, the input script should have a variable called 'result', which is the sole output of the
 * scriptlet execution.
 */
public class HostImpl extends Converter implements Host {
	private static boolean initPyEnvironment = false;
	
	static {
		if (!initPyEnvironment) {
			PySystemState.initialize();
		}
	}
	
	private PythonInterpreter interp = null;
	private String name = null;
	private ScriptLetParser p = null;
	private TreeMaker tm = new BaseTreeMaker();
	
	/**
	 * Constructor
	 */
	public HostImpl() throws IOException {
		name = "Python";
		
		p = new ScriptLetParser();
		p.setLexer(this);
		
		VersionImpl source_version = new VersionImpl(0, 7, 3);
		source = new LanguageVersionImpl("PharmML", source_version);
		
		VersionImpl target_version = new VersionImpl(2, 7, 6);
		target = new LanguageVersionImpl(name, target_version);
		
		interp = new PythonInterpreter();
	}
	
	@Override
	public BinaryTree createTree(Object o) { return tm.newInstance(o); }

	@Override
	public PyObject execute(String stmt) {
		if (stmt == null) throw new NullPointerException("The python syntax cannot be NULL.");
		interp.exec(stmt.toString()); 
		return interp.get("result");
	}
	
	@Override
	public IParser getParser() { return p; }

	@Override
	public PyObject getResult(String variable_name) {
		if (variable_name == null) return null;
		else return interp.get(variable_name);
	}

	
	@Override
	public String parse(Object context, Object element) {
		BinaryTree bt = tm.newInstance(element);
		for (NestedTreeRef nref : tm.getNestedTrees()) addStatement(nref);
		
		return p.parse(context, bt);
	}

	@Override
	public String getFunctionImpl(PyFunctionName ref) { return p.createPyFunction(ref); }
}

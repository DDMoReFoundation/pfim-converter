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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import crx.converter.spi.ILexer;
import crx.converter.spi.IParser;
import crx.converter.spi.ISettingReader;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.modellingsteps.OperationProperty;

/**
 * PFIM Settings Reader
 */
public class SettingReader implements ISettingReader {
	private Object ctx = new Object();
	private ILexer c = null;
	private IParser p = null;
	private Map<String, OperationProperty> prop_map = new HashMap<String, OperationProperty>();
	private boolean read_settings = false;
	List<OperationProperty> props = null; 
	
	private void addProperty(String name, OperationProperty prop) {
		if (name == null || prop == null) return;
		if (!prop_map.containsKey(name)) prop_map.put(name, prop);
	}
	
	@Override
	public Map<String, OperationProperty> getPropertyMap() { return prop_map; }

	String getValue(SettingLabel label) {
		if (p == null) return null;
		if (label == null) return null;
		if (hasValue(label)) {
			String value = p.parse(ctx, c.getStatement(prop_map.get(label.toString())));
			value = value.replaceAll("'", "");
			return value.trim();
		}
		
		return null;
	}
	
	@Override
	public String getValue(String name) {
		if (p == null) return null;
		if (hasValue(name)) {
			String value = p.parse(ctx, c.getStatement(prop_map.get(name)));
			value = value.replaceAll("'", "");
			return value.trim();
		}
		
		return null;
	}
 
	@Override
	public boolean hasReadSettings() { return read_settings; }

	boolean hasValue(SettingLabel label) {
		if (label == null) return false;
		else return hasValue(label.toString());
	}
	
	@Override
	public boolean hasValue(String name) {
		if (name == null) return false;
		else return prop_map.containsKey(name);
	}

	@Override
	public void readSettings() {
		if (c == null) return;
		if (props == null) return;	
		if (props.isEmpty()) return;
		
		TreeMaker tm = c.getTreeMaker();
		for (OperationProperty prop : props) {
			if (prop == null) continue;
			String name = prop.getName();
			if (name == null) continue;
			if (name.isEmpty()) continue;
				
			c.addStatement(prop, tm.newInstance(prop));
			c.updateNestedTrees();
			addProperty(name, prop);
		}
		
		read_settings = true;
	}

	@Override
	public void setLexer(ILexer c) { this.c = c; }

	@Override
	public void setParser(IParser p) { this.p = p; }

	
	@Override
	public void setProperties(List<OperationProperty> props) { this.props = props; }
}

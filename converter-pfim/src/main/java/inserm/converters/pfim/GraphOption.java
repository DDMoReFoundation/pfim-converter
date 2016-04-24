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

/**
 * Generic FIM calculation options.
 */
public enum GraphOption {
	GRAPH_LOGICAL("graphLogical"),
	GRAPH_SENSITIVITY("graphSensitivity"),
	LOG_GRAPH_OPTION("logarithmicGraphOption"),
	SAMPLING_TIME_LOWER("lowerSamplingTime"),
	SAMPLING_TIME_UPPER("upperSamplingTime"),
	STANDARD_GRAPHIC("standardGraphic"),
	X_AXES_NAMES("xAxesNames"),
	Y_AXES_NAMES("yAxesNames"),
	Y_AXES_RANGE("yAxesRange");
	
	public static boolean contains(String value){
		for (GraphOption item : values()) 
			if(item.toString().equals(value)) return true;
			
		return false;
	}
	
	private String value;
	
	private GraphOption(String value_){ value = value_; }
	
	public GraphOption fromValue(String value){
		for(GraphOption item : values())
			if (item.toString().equals(value)) return item;
			
		throw new IllegalArgumentException("Unknown enum type \""+ value+ "\".");
	}
	
	@Override
	public String toString(){ return value; }
}

/**
 * LOG graph options. 
 */
enum LogarithmicGraphOption {
	LOG_X("x"),
	LOG_XY("xy"),
	LOG_Y("y");
	
	private String value;
	
	private LogarithmicGraphOption(String value_){ value = value_; }
	
	public LogarithmicGraphOption fromValue(String value){
		for(LogarithmicGraphOption item : values())
			if (item.toString().equals(value)) return item;
			
		throw new IllegalArgumentException("Unknown enum type \""+ value+ "\".");
	}
	
	@Override
	public String toString(){ return value; }
}
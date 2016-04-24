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
 * Supported option for the Simplex algorithm.
 */
public enum SimplexOption {
	DELTA_TIME("deltaTime"),
	MAX_ITERATIONS("maximumIterations"),
	OPTIMISATION_OF_PROPORTIONS_OF_SUBJECTS("optimisationOfProportionsOfSubjects"),
	PRINT_ITERATIONS("printIterations"),
	RELATIVE_CONVERGENCE_TOLERANCE("relativeConvergenceTolerance"),
	SAMPLING_TIME_LOWER_A("lowerA"),
	SAMPLING_TIME_LOWER_B("lowerB"),
	SAMPLING_TIME_UPPER_A("upperA"),
	SAMPLING_TIME_UPPER_B("upperB"),
	SIMPLEX_PARAMETER("simplexParameter");
	
	public static boolean contains(String value) {
		for (SimplexOption item : values()) if (item.toString().equals(value)) return true;
		return false;
	}
	
	public static SimplexOption fromValue(String value){
		for(SimplexOption item : values())
			if (item.toString().equals(value)) return item;
			
		throw new IllegalArgumentException("Unknown enum type \""+ value+ "\".");
	}
	
	private String value;
	
	private SimplexOption(String value_){ value = value_; }
	
	@Override
	public String toString(){ return value; }
}

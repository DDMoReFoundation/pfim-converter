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
 * Supported PFIM error models.
 */
public enum SupportedErrorModel {
	PROPORTIONAL("proportionalErrorModel");
	
	public static boolean contains(String value) {
		for (SupportedErrorModel item : values()) if (item.toString().equals(value)) return true;
		return false;
	}
	
	private String value;
	
	private SupportedErrorModel(String value_){ value = value_; }
	
	public static SupportedErrorModel fromValue(String value){
		for(SupportedErrorModel item : values())
			if (item.toString().equals(value)) return item;
			
		throw new IllegalArgumentException("Unknown enum type \""+ value+ "\".");
	}
	
	@Override
	public String toString(){ return value; }
}

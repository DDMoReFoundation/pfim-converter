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

import java.io.File;
import java.io.IOException;

import crx.converter.engine.Lexer;
import crx.converter.engine.parts.EstimationStepImpl;
import eu.ddmore.convertertoolbox.domain.LanguageVersionImpl;
import eu.ddmore.convertertoolbox.domain.VersionImpl;

public class Converter extends Lexer {
	/**
	 * Major version
	 */
	public static final String MAJOR_VERSION = "4";
	
	/**
	 * Converter name
	 */
	public static final String NAME = "PFIM";
	
	public Converter() throws IOException, NullPointerException {
		super();
		name = NAME;
		
		// Register lexer/parser dependency.
		Parser p = new Parser();
		setParser(p);
		p.setLexer(this);
		
		source = createPharmMLVersion();
		
		VersionImpl target_version = new VersionImpl(Integer.parseInt(MAJOR_VERSION), 0, 0);
		target = new LanguageVersionImpl("PFIM", target_version);
	}
	
	@Override
	//  Overriding this method in case Java generates Windows style file paths that 'R' does not like.
	public String getOutputDirectory() {
		if (output_dir == null) return null;
		
		try {
			File dir = new File(output_dir);
			output_dir = dir.getCanonicalPath();
			output_dir = output_dir.replace('\\', '/');
		} catch (Exception e) {
			if (is_echo_exception) e.printStackTrace();
		}
		
		return output_dir; 
	}
	
	@Override
	protected void initialise() {
		EstimationStepImpl.setUseDefaultParameterEstimate(true);
		EstimationStepImpl.setDefaultParameterEstimateValue(1.0);
		setUsePiecewiseAsEvents(true);
		sort_parameter_model = true;
		sort_structural_model = true;
	}
}

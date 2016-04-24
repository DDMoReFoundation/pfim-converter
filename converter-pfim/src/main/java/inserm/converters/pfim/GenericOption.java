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
 * Generic option flags for PFIM.
 */
public enum GenericOption {
	ALPHA("alpha"),
	ATOL("atol"),
	BETA_COVARIATE("betaCovariate"),
	BETA_COVARIATE_OCCASSION("occassionBetaCovariate"),
	COMPUTE_NNI("computeNNI"),
	COMPUTE_NNI_EQUIVALENCE("computeNNIEquivalence"),
	COMPUTE_POWER("computePower"),
	COMPUTE_POWER_EQUIVALENCE("computePowerEquivalence"),
	COVARIATE_OCCASSION_CATEGORIES("covariateOccassionCategories"),
	COVARIATE_OCCASSION_NAMES("covariateOccassionNames"),
	COVARIATE_OCCASSION_PROPORTIONS("covariateOccassionProportions"),
	COVARIATE_OCCASSION_SEQUENCE("covariateOccassionSequence"),
	COVARIATE_PROPORTIONS("covariateProportions"),
	DOSE_IDENTICAL("doseIdentical"),
	DOSE_VECTOR("doseVector"),
	EQUIVALENCE_INTERVAL_START("equivalenceIntervalStart"),
	EQUIVALENCE_INTERVAL_STOP("equivalenceIntervalStop"),
	FULL_FIM("calculateFullFIM"),
	GAMMA("gamma"),
	GIVEN_POWER("givenPower"),
	IDENTICAL_COND_4_EACH_DESIGN("condinitIdentical"),
	IDENTICAL_TIMES("identicalTimes"),
	N_OCC("numberOfOccassions"),
	NR("numberOfResponses"),
	NSUBJECTS("numberOfSubjects"),
	OMEGA("omega"),
	PARAMETER_ASSOCIATED("parameterAssociated"),
	PARAMS_ASSOCIATED_WITH_OCCASSION("paramsAssociatedWithOccassion"),
	PROPORTIONS("usingProportions"),
	PROTA("protA"),
	RTOL("rtol"),
	SIG_INTERA("sig.interA"),
    SIG_SLOPEA("sig.slopeA"),
    TRAND("randomEffectModel"),
    USING_COVARIATE("usingCovariateModel"),
    USING_COVARIATE_OCCASSION("useCovariateOccassionModel");
	
	public static boolean contains(String value) {
		if (value == null) return false;
		if (value.startsWith(COVARIATE_PROPORTIONS.toString())) return true;
		else if (value.startsWith(PARAMETER_ASSOCIATED.toString())) return true;
		else if (value.startsWith(BETA_COVARIATE.toString())) return true;
		else if (value.startsWith(COVARIATE_OCCASSION_CATEGORIES.toString())) return true;
		else if (value.startsWith(COVARIATE_OCCASSION_SEQUENCE.toString())) return true;
		else if (value.startsWith(COVARIATE_OCCASSION_PROPORTIONS.toString())) return true;
		else if (value.startsWith(PARAMS_ASSOCIATED_WITH_OCCASSION.toString())) return true;
		else if (value.startsWith(BETA_COVARIATE_OCCASSION.toString())) return true;
		
		for (GenericOption item : values()) 
			if (item.toString().equals(value)) return true;
			
		return false;
	}
	
	public static GenericOption fromValue(String value){
		for(GenericOption item : values())
			if (item.toString().equals(value)) return item;
			
		throw new IllegalArgumentException("Unknown enum type \""+ value+ "\".");
	}
	
	private String value;
	
	private GenericOption(String value_){ value = value_; }
	
	@Override
	public String toString(){ return value; }
}

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

package inserm.converters.pfim.parts;

import java.util.ArrayList;
import java.util.List;

import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.modeldefn.CovariateDefinition;
import eu.ddmore.libpharmml.dom.modeldefn.IndividualParameter;
import eu.ddmore.libpharmml.dom.modeldefn.TransformationType;

/**
 * A fixed effect reference in a covariate referenced in PFIM. 
 */
public class FixedEffectRef extends PharmMLRootType {
	public List<PharmMLRootType> beta_covariates = new ArrayList<PharmMLRootType>();
	public CovariateDefinition cov = null;
	public List<IndividualParameter> ips = new ArrayList<IndividualParameter>();
	public TransformationType transformation = TransformationType.IDENTITY;
	
	/**
	 * Constructor
	 * @param cov_ Covariate
	 */
	public FixedEffectRef(CovariateDefinition cov_) {
		if (cov_ == null) throw new NullPointerException("Covariate cannot be NULL");
		cov = cov_;
	}
	
	/**
	 * Add a beta covariate to the fixed effect covariate reference.
	 * This is normally a numeric value linked to a population parameter.
	 * @param element Parameter
	 * @return boolean
	 */
	public boolean addBetaCovariate(PharmMLRootType element) {
		if (element != null) {
			if (!beta_covariates.contains(element)) {
				beta_covariates.add(element);
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Add an individual parameter to the fixed effect covariate reference.
	 * @param ip Individual Parameter
	 * @return boolean
	 */
	public boolean addIndividualParameter(IndividualParameter ip) {
		if (ip != null) {
			if (!ips.contains(ip)) {
				ips.add(ip);
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Check if the fixed effect reference has Beta Covariates.
	 * @return boolean
	 */
	public boolean hasBetaCovariates() { return beta_covariates.size() > 0; }
	
	/**
	 * Check if the fixed effect reference has individual parameters.
	 * @return boolean
	 */
	public boolean hasIndividualParameter() { return ips.size() > 0; }
	
	@Override 
	public String toString() {
		String format = "%s=%s";
		StringBuffer sb = new StringBuffer("FixedEffectRef:{");
		if (cov != null) {
			sb.append(String.format(format, "covariate", cov.getSymbId()));
			sb.append(",");
		}
		sb.append(String.format(format, "transformation", transformation));
		sb.append(",");
		
		List<String> params = new ArrayList<String>();
		for (IndividualParameter ip : ips) {
			if (ip == null) params.add(null);
			else params.add(ip.getSymbId());
		}
		
		sb.append(String.format(format, "parameters", params));
		sb.append(",");
		
		sb.append(String.format(format, "beta_covariates", beta_covariates));
		sb.append("}");
		return sb.toString();
	}
}

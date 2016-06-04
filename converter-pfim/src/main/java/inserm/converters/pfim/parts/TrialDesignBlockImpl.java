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

import static crx.converter.engine.PharmMLTypeChecker.isBolus;
import static crx.converter.engine.PharmMLTypeChecker.isObservation;
import static crx.converter.engine.Utils.getClassName;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.Accessor;
import crx.converter.engine.common.DataFiles;
import crx.converter.spi.ILexer;
import crx.converter.spi.IParser;
import crx.converter.spi.blocks.TrialDesignBlock;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.OidRef;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.StandardAssignable;
import eu.ddmore.libpharmml.dom.trialdesign.Administration;
import eu.ddmore.libpharmml.dom.trialdesign.ArmDefinition;
import eu.ddmore.libpharmml.dom.trialdesign.Arms;
import eu.ddmore.libpharmml.dom.trialdesign.Bolus;
import eu.ddmore.libpharmml.dom.trialdesign.DosingRegimen;
import eu.ddmore.libpharmml.dom.trialdesign.DosingVariable;
import eu.ddmore.libpharmml.dom.trialdesign.InterventionSequence;
import eu.ddmore.libpharmml.dom.trialdesign.Interventions;
import eu.ddmore.libpharmml.dom.trialdesign.Observation;
import eu.ddmore.libpharmml.dom.trialdesign.ObservationSequence;
import eu.ddmore.libpharmml.dom.trialdesign.Observations;
import eu.ddmore.libpharmml.dom.trialdesign.TrialDesign;

/**
 * Wrapper class for the PharmML trial design block.
 */
public class TrialDesignBlockImpl extends PartImpl implements TrialDesignBlock {
	public static class InterventionSequenceRef {
		public String administration_oid = null;
		public double start = 0.0;
		
		public InterventionSequenceRef(OidRef admin_ref, double start_) {
			if (admin_ref == null) throw new NullPointerException("OID Reference is NULL");
			administration_oid = admin_ref.getOidRef();
			if (administration_oid == null) throw new NullPointerException("OID is NULL");
			start = start_;
		}
	}
	
	private Map<ArmDefinition,String> arm_2_observation_map = new  HashMap<ArmDefinition,String>(); 
	private Map<String, ArmDefinition> arm_map = new HashMap<String, ArmDefinition>();
	private Map<ArmDefinition,Double> arm_observation_start_map = new  HashMap<ArmDefinition,Double>();
	private Map<String, Integer> arm_size_map = new HashMap<String, Integer>();
	private List<ArmDefinition> arms = new ArrayList<ArmDefinition>();
	private Object ctx = new Object();
	private Map<String, String> dose_stmt_map = new HashMap<String, String>();
	private Map<String, PharmMLRootType> dose_target_map = new HashMap<String, PharmMLRootType>();
	private Map<String, Double> dose_time_raw = new HashMap<String, Double>();
	private Map<String, InterventionSequenceRef> iseq_map = new HashMap<String, InterventionSequenceRef>(); 
	private List<Observation> obs = new ArrayList<Observation>();
	private Map<String, Observation> obs_map = new HashMap<String, Observation>();
	private IParser p = null;
	private TrialDesign td = null;
	private TreeMaker tm = null;
	
	/**
	 * Constructor
	 * @param td_ Trial Design Block
	 * @param c_ Converter Instance
	 */
	public TrialDesignBlockImpl(TrialDesign td_, ILexer c_) {
		if (td_ == null) throw new NullPointerException("The trial design object is NULL.");
		if (c_ == null) throw new NullPointerException("The converter object is NULL.");
		
		td = td_;
		lexer = c_;
		tm = c_.getTreeMaker();
		p = c_.getParser();
		
		DataFiles dfs = lexer.getDataFiles();
		if (dfs != null) dfs.setExternalDataSets(td.getListOfExternalDataSet()); 
	}
	
	private void addArm(ArmDefinition arm) {
		String oid = arm.getOid();
		if (oid == null) throw new NullPointerException("Arm identifier is NULL");
		if (!arm_map.containsKey(oid)) arm_map.put(oid, arm);
	} 
	private void buildArms() {
		Arms arm_list = td.getArms();
		if (arm_list == null) return;
		
		for (ArmDefinition arm : arm_list.getListOfArm()) {
			if (arm == null)  continue;
			arms.add(arm);
			addArm(arm);
			processSize(arm);
			processInterventionSequence(arm);
			processObservationSequence(arm);
		}
	}
	
	private void buildInterventions() {
		Interventions int_list = td.getInterventions();
		if (int_list == null) return;
		
		for (Administration admin : int_list.getListOfAdministration()) {
			if (admin == null) continue;
			String oid = admin.getOid();
			if (oid == null) throw new NullPointerException("Administration OID is not set");
			JAXBElement<? extends DosingRegimen> tag = admin.getDosingRegimen();
			if (tag == null) continue;
			DosingRegimen regimen = tag.getValue();
			if (isBolus(regimen)) processBolus(oid, (Bolus) regimen);
			else throw new UnsupportedOperationException("Unsupported dosing regimen (class='" + getClassName(regimen) + "'");
		}
	} 
	
	private void buildObservations() {
		Observations observations = td.getObservations();
		if (observations == null) return;
		for (PharmMLRootType element : observations.getListOfObservationsElements()) {
			if (isObservation(element)) {
				Observation ob = (Observation) element;
				String oid = ob.getOid();
				if (oid == null) throw new NullPointerException("The observation OID is NULL");
				if (obs_map.containsKey(oid)) throw new IllegalStateException("An observation OID is not unique (oid='" + oid + "')");
				obs.add(ob);
				obs_map.put(oid, ob);
				
				lexer.addStatement(ob.getObservationTimes(), tm.newInstance(ob.getObservationTimes()));
				lexer.updateNestedTrees();
			}
		}
	}
	
	@Override
	public void buildTrees() {
		buildObservations();
		buildInterventions();
		buildArms();
	}
	
	/**
	 * Get the dose time associated with an ad
	 * @param admin_oid
	 * @return
	 */
	public double getAdministrationStartTime(String admin_oid) {
		if (admin_oid == null) return 0.0;
		if (dose_time_raw.containsKey(admin_oid)) return dose_time_raw.get(admin_oid);
		else return 0.0;
	}
	
	@Override
	public int getArmCount() { return arms.size(); }
	
	@Override
	public List<ArmDefinition> getArms() { return arms; };
	
	@Override
	public int getArmSize(String oid) {
		if (oid == null) return 0;
		int size = 0;
		if (arm_size_map.containsKey(oid)) size = arm_size_map.get(oid);
		return size;
	}
	
	@Override
	public String getDoseStatement(String administration_oid) {
		if (administration_oid == null) return null;
		if (dose_stmt_map.containsKey(administration_oid)) return dose_stmt_map.get(administration_oid);
		else return null;
	} 
	
	@Override
	public PharmMLRootType getDoseTarget(String administration_oid) {
		if (administration_oid == null) return null;
		if (dose_target_map.containsKey(administration_oid)) return dose_target_map.get(administration_oid);
		else return null;
	}
	
	/**
	 * Get the intervention reference linked to the arm.
	 * @param arm Arm Definition
	 * @return InterventionSequenceRef
	 */
	public InterventionSequenceRef getInterventionSequenceRef(ArmDefinition arm) {
		if (arm == null) return null;
		String oid = arm.getOid();
		if (iseq_map.containsKey(oid)) return iseq_map.get(oid);
		else return null;
	}
	
	/**
	 * Get the model/source for the trial design block.
	 * @return TrialDesign
	 */
	public TrialDesign getModel() { return td; }
	
	@Override
	public String getName() { return "trial_design"; }
	
	/**
	 * Get the observation linked to the Arm.
	 * @param arm Arm Instance
	 * @return Observation
	 */
	public Observation getObservation(ArmDefinition arm) {
		if (arm == null) return null;
		if (arm_2_observation_map.containsKey(arm)) {
			String ob_oid = arm_2_observation_map.get(arm);
			return this.getObservation(ob_oid);
		}
		
		return null;
	}
	
	/**
	 * Get the named observation element
	 * @param oid Observation Identifier
	 * @return Observation
	 */
	public Observation getObservation(String oid) {
		if (oid == null) return null;
		else if (obs_map.containsKey(oid))  return obs_map.get(oid);
		else return null;
	}
	
	/**
	 * Get list of declared observations.
	 * @return
	 */
	public List<Observation> getObservations() { return obs; }
	
	/**
	 * Get the start offset for an Arm
	 * @param arm Arm Instance
	 * @return double
	 */
	public double getObservationStart(ArmDefinition arm) {
		double start = 0.0;
		if (arm != null) {
			if (arm_observation_start_map.containsKey(arm)) start = arm_observation_start_map.get(arm);
		}
		
		return start;
	}
	
	@Override
	public List<DerivativeVariable> getStateVariablesWithDosing() { throw new UnsupportedOperationException(); }

	@Override
	public List<String> getSymbolIds() { return new ArrayList<String>(); }

	/**
	 * Flag if a state variable is associated with a dosing event.
	 * @return boolean
	 */
	public boolean hasDosing() {
		return getStateVariablesWithDosing().isEmpty() == false;
	}
	
	@Override
	public boolean hasOccassions() { return td.getListOfOccasions().size() > 0; }
	
	@Override
	public boolean hasSymbolId(String name) { return false; }
	
	private void processBolus(String oid, Bolus bolus) {
		if (oid == null || bolus == null) return;
		
		DosingVariable target = bolus.getDoseAmount();
		if (target == null) throw new NullPointerException("Dose target not specified (oid='" + oid + "')");
		Accessor a = lexer.getAccessor();
		
		PharmMLRootType element = a.fetchElement(target);
		if (element == null) throw new NullPointerException("Dose target element not found in the model (oid='" + oid + "')");
		dose_target_map.put(oid, element);
		
		// Assuming that the dose amount is an numeric quantity.
		String stmt = stripOuterBrackets(p.parse(ctx, tm.newInstance(target)).trim());
		dose_stmt_map.put(oid, stmt);
		BinaryTree bt = tm.newInstance(bolus.getDosingTimes());
		lexer.updateNestedTrees();
		dose_time_raw.put(oid, Double.parseDouble(p.parse(ctx, bt).trim()));
	}
	
	private void processInterventionSequence(ArmDefinition arm) {
		if (arm == null) return;
		String oid = arm.getOid();
		if (oid == null) throw new NullPointerException("OID is undefined");
		
		List<InterventionSequence> iseqs = arm.getListOfInterventionSequence();
		if (iseqs.size() == 0) throw new IllegalStateException("No intervention sequence specified.");
		if (iseqs.size() > 1) throw new IllegalStateException("Multiple intervention sequences specified.");
		
		double start = 0.0;
		InterventionSequence iseq = iseqs.get(0);
		if (iseq.getStart() != null) start = Double.parseDouble(p.parse(ctx, tm.newInstance(iseq.getStart())).trim());
		
		OidRef oref = iseq.getInterventionList().getListOfInterventionRef().get(0);
		iseq_map.put(oid, new InterventionSequenceRef(oref, start));
	}
	
	// As PFIM, assuming 1 observation sequence per arm.
	private void processObservationSequence(ArmDefinition arm) {
		if (arm == null) return;
		ObservationSequence oseq = arm.getListOfObservationSequence().get(0);
		if (oseq == null) return;
		OidRef oref = oseq.getObservationList().getListOfObservationRef().get(0);
		if (oref == null) return;
		if (oref.getOidRef() == null) return;
		arm_2_observation_map.put(arm, oref.getOidRef());
		
		Double start = 0.0;
		if (oseq.getStart() != null) 
			start = Double.parseDouble(p.parse(ctx, tm.newInstance(oseq.getStart())).trim());
		arm_observation_start_map.put(arm, start);
	}
	
	private void processSize(ArmDefinition arm) {
		if (arm == null) return;
		if (arm.getOid() == null) return;
		
		StandardAssignable size_expr = arm.getArmSize();
		if (size_expr == null) return;
		
		p = lexer.getParser();
		Integer size  = Integer.parseInt(p.parse(ctx, tm.newInstance(size_expr)).trim());
		arm_size_map.put(arm.getOid(), size);
	}
	
	private String stripOuterBrackets(String stmt) {
		if (stmt == null) return null;
		stmt = stmt.trim();
		
		if (stmt.startsWith("(") && stmt.endsWith(")")) {
			int indexOfOpenBracket = stmt.indexOf("(");
			int indexOfLastBracket = stmt.lastIndexOf(")");
			stmt = stmt.substring(indexOfOpenBracket+1, indexOfLastBracket);
		}
		
		return stmt;
	}
}

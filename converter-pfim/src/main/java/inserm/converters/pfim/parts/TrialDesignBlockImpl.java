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
import static crx.converter.engine.PharmMLTypeChecker.isStructuredError;
import static crx.converter.engine.Utils.getClassName;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBElement;

import crx.converter.engine.Accessor;
import crx.converter.engine.common.DataFiles;
import crx.converter.engine.common.ElementaryDesign;
import crx.converter.engine.common.InterventionSequenceRef;
import crx.converter.engine.common.Protocol;
import crx.converter.spi.ILexer;
import crx.converter.spi.IParser;
import crx.converter.spi.blocks.ObservationBlock;
import crx.converter.spi.blocks.TrialDesignBlock2;
import crx.converter.tree.BinaryTree;
import crx.converter.tree.TreeMaker;
import eu.ddmore.libpharmml.dom.commontypes.DerivativeVariable;
import eu.ddmore.libpharmml.dom.commontypes.OidRef;
import eu.ddmore.libpharmml.dom.commontypes.PharmMLRootType;
import eu.ddmore.libpharmml.dom.commontypes.StandardAssignable;
import eu.ddmore.libpharmml.dom.modeldefn.StructuredObsError;
import eu.ddmore.libpharmml.dom.trialdesign.Administration;
import eu.ddmore.libpharmml.dom.trialdesign.ArmDefinition;
import eu.ddmore.libpharmml.dom.trialdesign.Arms;
import eu.ddmore.libpharmml.dom.trialdesign.Bolus;
import eu.ddmore.libpharmml.dom.trialdesign.DesignSpaces;
import eu.ddmore.libpharmml.dom.trialdesign.DosingRegimen;
import eu.ddmore.libpharmml.dom.trialdesign.DosingVariable;
import eu.ddmore.libpharmml.dom.trialdesign.InterventionSequence;
import eu.ddmore.libpharmml.dom.trialdesign.Interventions;
import eu.ddmore.libpharmml.dom.trialdesign.Observation;
import eu.ddmore.libpharmml.dom.trialdesign.ObservationList;
import eu.ddmore.libpharmml.dom.trialdesign.ObservationSequence;
import eu.ddmore.libpharmml.dom.trialdesign.Observations;
import eu.ddmore.libpharmml.dom.trialdesign.ObservationsCombination;
import eu.ddmore.libpharmml.dom.trialdesign.SingleDesignSpace;
import eu.ddmore.libpharmml.dom.trialdesign.SingleObservation;
import eu.ddmore.libpharmml.dom.trialdesign.TrialDesign;

/**
 * Wrapper class for the PharmML trial design block.
 */
public class TrialDesignBlockImpl extends PartImpl implements TrialDesignBlock2 {
	private Accessor a = null; 
	private Map<String, Administration> admin_map_ = new HashMap<String, Administration>();
	private List<Administration> admins = new ArrayList<Administration>();
	private Map<ArmDefinition,String> arm_2_observation_map = new  HashMap<ArmDefinition,String>();
	private Map<String, ArmDefinition> arm_map = new HashMap<String, ArmDefinition>();
	private Map<Observation, List<ArmDefinition>> arm_membership_map = new HashMap<Observation, List<ArmDefinition>>();
	private Map<ArmDefinition,Double> arm_observation_start_map = new  HashMap<ArmDefinition,Double>();
	private Map<String, Integer> arm_size_map = new HashMap<String, Integer>();
	private List<ArmDefinition> arms = new ArrayList<ArmDefinition>();
	private Object ctx = new Object(); 
	private List<SingleDesignSpace> design_spaces = new ArrayList<SingleDesignSpace>();
	private Map<SingleDesignSpace, Observation> designspace_2_window = new HashMap<SingleDesignSpace, Observation>();
	private Map<String, String> dose_stmt_map = new HashMap<String, String>();
	private Map<String, PharmMLRootType> dose_target_map = new HashMap<String, PharmMLRootType>();
	private Map<String, Double> dose_time_raw = new HashMap<String, Double>();
	private Map<String, InterventionSequenceRef> iseq_map = new HashMap<String, InterventionSequenceRef>();
	private List<Observation> obs = new ArrayList<Observation>();
	private Map<String, Observation> obs_map = new HashMap<String, Observation>();
	private IParser p = null; 
	private List<Protocol> protocols = new ArrayList<Protocol>();  
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
		a = c_.getAccessor();
		
		DataFiles dfs = lexer.getDataFiles();
		if (dfs != null) dfs.setExternalDataSets(td.getListOfExternalDataSet()); 
	} 
	
	private void addArm(ArmDefinition arm) {
		String oid = arm.getOid();
		if (oid == null) throw new NullPointerException("Arm identifier is NULL");
		if (!arm_map.containsKey(oid)) arm_map.put(oid, arm);
	}
	
	private void addArmToList(ArmDefinition arm, List<ArmDefinition> list) {
		if (arm == null || list == null) return;
		if (!list.contains(arm)) list.add(arm);
	} 
	
	private void buildAdministrations_() {
		Interventions ints = td.getInterventions();
		if (ints == null) return;
		
		List<Administration> admins_ = ints.getListOfAdministration();
		if (admins_ == null) return;
		if (admins_.isEmpty()) return;
	
		for (Administration admin : admins_) {
			if (admin == null) continue;
			String oid = admin.getOid();
			if (oid == null) throw new NullPointerException("Administration OID is NULL");
			admins.add(admin);
			admin_map_.put(oid, admin);
		}
	}
	
	private void buildArmMembershipMap() {
		initArmMembershipMap();
		
		Map<String, ObservationsCombination> comb_map = createObservationsCombinationMap();
		Map<ObservationsCombination, List<Observation>> comb_2_obs_map = createObservationCombinbination2ObservationMap();
		
		for (ArmDefinition arm : arms) {
			if (arm == null) continue;
			List<ObservationSequence> oseqs = arm.getListOfObservationSequence();
			if (oseqs == null) continue;
			for (ObservationSequence oseq : oseqs) {
				if (oseq == null) continue;
				ObservationList outer_list = oseq.getObservationList();
				if (outer_list == null) continue;
				List<OidRef> ob_refs = outer_list.getListOfObservationRef();
				if (ob_refs == null) continue;
				for (OidRef ob_ref : ob_refs) {
					if (ob_ref == null) continue;
					String oid = ob_ref.getOidRef();
					if (oid == null) continue;
					
					// Look in Observation combination
					if (comb_map.containsKey(oid)) {
						ObservationsCombination comb = comb_map.get(oid);
						if (comb != null) {
							List<Observation> referenced_obs = comb_2_obs_map.get(comb);
							if (referenced_obs != null) {
								for (Observation referenced_ob : referenced_obs) {
									if (referenced_ob == null) continue;
									if (arm_membership_map.containsKey(referenced_ob)) 
										addArmToList(arm, arm_membership_map.get(referenced_ob));
								}
							}
						}
					} else {
						Observation referenced_ob = getObservation(oid);
						if (referenced_ob != null) {
							if (arm_membership_map.containsKey(referenced_ob)) 
								addArmToList(arm, arm_membership_map.get(referenced_ob));
						}
					}
				}
			}
		}
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
	
	private void buildDesignSpaces() {
		DesignSpaces ds = td.getDesignSpaces();
		if (ds == null) return;
		
		processDesignSpaces(ds.getListOfDesignSpace());
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
	
	private void buildProtocols() {
		List<ObservationBlock> oblocks = lexer.getObservationBlocks();
		if (oblocks.isEmpty()) return;
		
		Map<ObservationBlock, PharmMLRootType> output_map = createOutputMap();
		Map<ObservationBlock, PharmMLRootType> standard_map = createStandardMap();
		
		List<String> windows = new ArrayList<String>();
		Map<String, List<PharmMLRootType>> obs_2_element_map = createObservationToElementMap(windows);
		Map<String, List<StandardAssignable>> ob_timepoint_map = createObservationToTimepointMap(windows, obs_2_element_map);
		
		// Generate the protocol specification for each response.
		int i = 0;
		for (ObservationBlock oblock : oblocks) {
			Protocol protocol = createProtocol(oblock, output_map, obs_2_element_map, ob_timepoint_map, windows); 
			if (protocol == null) protocol = createProtocol(oblock, standard_map, obs_2_element_map, ob_timepoint_map, windows);
			
			if (protocol != null) {
				protocol.index = i++;
				protocols.add(protocol);
			} else
				throw new IllegalStateException("Unable to derive a protocol for observation model (oid='" + oblock + "')");
		}
	}
	
	private void buildStartTimes() {
		for (Protocol protocol : protocols) {
			if (protocol == null) continue;
			for (ElementaryDesign ed : protocol.elementary_designs) {
				if (ed == null) continue;
				if (ed.observation_oid == null) continue;
				
				ArmDefinition arm = getArm(ed.arm_oid);
				if (arm == null) continue;
				if (arm.getListOfObservationSequence() != null) {
					for (ObservationSequence oseq : arm.getListOfObservationSequence()) {
						if (oseq == null) continue;
				
						ObservationList referenced_windows = oseq.getObservationList();
						if (referenced_windows == null) continue;
						if (referenced_windows.getListOfObservationRef() == null) continue;
						
						boolean foundObservation = false;
						for (OidRef ref : referenced_windows.getListOfObservationRef()) {
							if (ref == null) continue;
							String oid = ref.getOidRef();
							if (ed.observation_oid.equals(oid)) {
								foundObservation = true;
								break;
							}
						}
						
						// Register a start time with the elementary design.
						if (oseq.getStart() != null && foundObservation) 
							ed.start_time_offset = 
								Double.parseDouble(p.parse(ctx, tm.newInstance(oseq.getStart())).trim());
					}
				}
			}
		}
	}
	
	@Override
	public void buildTrees() {
		buildAdministrations_();
		buildObservations();
		buildInterventions();
		buildArms();
		buildDesignSpaces();
		buildArmMembershipMap();
		buildProtocols();
		buildStartTimes();
	}
	
	private Map<ObservationsCombination, List<Observation>> createObservationCombinbination2ObservationMap() {
		Map<ObservationsCombination, List<Observation>> map = new HashMap<ObservationsCombination, List<Observation>>();
		
		Observations observations = td.getObservations();
		if (observations == null) return map;
		List<ObservationsCombination> ob_combs = observations.getListOfObservationsCombination();
		if (ob_combs == null) return null;
		for (ObservationsCombination ob_comb : ob_combs) {
			List<Observation> referenced_obs = new ArrayList<Observation>();
			if (ob_comb == null) continue;
			List<SingleObservation> single_obs = ob_comb.getListOfObservations();
			if (single_obs != null) {
				for (SingleObservation single_ob : single_obs) {
					if (single_ob == null) continue;
					if (single_ob.getListOfObservationRef() != null) {
						for (OidRef ref : single_ob.getListOfObservationRef()) {
							if (ref == null) continue;
							Observation ob = getObservation(ref);
							if (ob != null ) if (!referenced_obs.contains(ob)) referenced_obs.add(ob);
						}
					}
				}
			}
			
			map.put(ob_comb, referenced_obs);
		}
		
		return map;
	}
	
	private Map<String, ObservationsCombination> createObservationsCombinationMap() {
		Map<String, ObservationsCombination> map = new HashMap<String, ObservationsCombination>();
		
		Observations observations = td.getObservations();
		if (observations == null) return map;
		List<ObservationsCombination> ob_combs = observations.getListOfObservationsCombination();
		if (ob_combs == null) return null;
		
		for (ObservationsCombination ob_comb : ob_combs) {
			if (ob_comb == null) continue;
			String oid = ob_comb.getOid();
			if (oid == null) throw new NullPointerException("An observation combination was not assigned an OID");
			if (map.containsKey(oid)) throw new IllegalStateException("An observation combination OID is not unique");
			map.put(oid, ob_comb);
		}
		
		return map;
	}
	
	Map<String, List<PharmMLRootType>> createObservationToElementMap(List<String> windows) {
		Map<String, List<PharmMLRootType>> obs_2_element_map = new HashMap<String, List<PharmMLRootType>>();
		if (windows != null) {
			for (Observation ob : obs) {
				String oid = ob.getOid();
				if (oid == null) throw new NullPointerException("An observation OID is not specified");
				List<PharmMLRootType> elements = a.fetchElementList(ob.getContinuous());
				obs_2_element_map.put(oid, elements);
				if (!windows.contains(oid)) windows.add(oid);
			}
		}
		
		return obs_2_element_map;
	}
	
	Map<String, List<StandardAssignable>> createObservationToTimepointMap(List<String> windows, Map<String, List<PharmMLRootType>> obs_2_element_map) {
		Map<String, List<StandardAssignable>> ob_timepoint_map = new HashMap<String, List<StandardAssignable>>();
		if (windows != null && obs_2_element_map != null) {
			Observations observations = td.getObservations();
			if (observations == null) return ob_timepoint_map;

			// List elements the are associated possibly with a combined sampling windows.
			if (observations.getListOfObservationsCombination() != null) {
				for (ObservationsCombination oc : observations.getListOfObservationsCombination()) {
					if (oc == null) continue;
					List<PharmMLRootType> elements = new ArrayList<PharmMLRootType>();
					if (oc.getListOfObservations() == null) continue;
					String oid = oc.getOid(); 
					if (oid == null) throw new NullPointerException("An observation combination identifier is NULL");
					if (!windows.contains(oid)) windows.add(oid);
					for (SingleObservation so : oc.getListOfObservations()) {
						if (so == null) continue;
						if (so.getListOfObservationRef() == null) continue;
						for (OidRef ref : so.getListOfObservationRef()) {
							if (ref == null) continue;
							String inner_oid = ref.getOidRef();
							if (inner_oid == null) continue;
							if (obs_2_element_map.containsKey(inner_oid)) {
								List<PharmMLRootType> other_elements = obs_2_element_map.get(inner_oid);
								if (other_elements != null) elements.addAll(other_elements);
							}
						}
					}

					obs_2_element_map.put(oid, elements);
					if (!windows.contains(oid)) windows.add(oid);
				}
			}

			// Create time point maps for the registered sampling windows,
			if (observations.getListOfObservationsElements() != null) {
				for (PharmMLRootType o : observations.getListOfObservationsElements()) {
					if (o == null) continue;
					if (isObservation(o)) {
						Observation ob = (Observation) o;
						String oid = ob.getOid();
						if (oid == null) throw new NullPointerException("An observation OID is NULL");
						if (ob.getObservationTimes() != null) {
							List<StandardAssignable> tps = new ArrayList<StandardAssignable>();
							tps.add(ob.getObservationTimes());
							ob_timepoint_map.put(oid, tps);
						}
					}
				}
			}

			// Register sampling expressions for each window.
			if (observations.getListOfObservationsCombination() != null) {
				for (ObservationsCombination oc : observations.getListOfObservationsCombination()) {
					if (oc == null) continue;
					if (oc.getListOfObservations() == null) continue;
					if (oc.getOid() == null) throw new NullPointerException("An observation combination identifier is NULL");
					List<StandardAssignable> tps = new ArrayList<StandardAssignable>();
					for (SingleObservation so : oc.getListOfObservations()) {
						if (so == null) continue;
						if (so.getListOfObservationRef() == null) continue;
						for (OidRef ref : so.getListOfObservationRef()) {
							if (ref == null) continue;
							String oid = ref.getOidRef();
							if (oid == null) continue;
							if (ob_timepoint_map.containsKey(oid)) {
								List<StandardAssignable> other_tps = ob_timepoint_map.get(oid);
								for (StandardAssignable other_tp : other_tps) {
									if (other_tp == null) continue;
									if (!tps.contains(other_tp)) tps.add(other_tp);
								}
							}
						}
					}
					ob_timepoint_map.put(oc.getOid(), tps);
				}
			}
		}

		return ob_timepoint_map;
	}
	
	private Map<ObservationBlock, PharmMLRootType> createOutputMap() {
		Map<ObservationBlock, PharmMLRootType> output_map = new HashMap<ObservationBlock, PharmMLRootType>(); 
		List<ObservationBlock> oblocks = lexer.getObservationBlocks();
		if (oblocks.isEmpty()) return output_map;
		
		for (ObservationBlock block : oblocks) {
			if (block == null) continue;
		
			if (!isStructuredError(block.getObservationError())) 
				throw new IllegalStateException("Unable to build a protocol as error model is not a structured error model.");
			
			StructuredObsError soe = (StructuredObsError) block.getObservationError();
			PharmMLRootType element = a.fetchElement(soe.getOutput());
			if (element == null) throw new NullPointerException("Model element referenced in an error model was not found");
			output_map.put(block, element);
			
		}
		
		return output_map;
	}
	
	private Protocol createProtocol(ObservationBlock oblock, 
			Map<ObservationBlock, PharmMLRootType> block_map, Map<String, 
			List<PharmMLRootType>> obs_2_element_map,
			Map<String, List<StandardAssignable>> ob_timepoint_map, List<String> windows) 
	{
		Protocol protocol = new Protocol(oblock, 0);
		PharmMLRootType referenced_output = block_map.get(oblock);
		if (referenced_output == null) 
			throw new NullPointerException("Protocol output variable not specified (oid='" + oblock + "')");

		for (String window : windows) {
			if (window == null) continue;
			List<PharmMLRootType> referenced_elements = obs_2_element_map.get(window);
			
			boolean hasAssociation = false;
			if (referenced_elements != null) {
				for (PharmMLRootType referenced_element : referenced_elements) {
					if (referenced_element == null) continue;
					if (referenced_output.equals(referenced_element)) {
						hasAssociation = true;
						break;
					}
				}
			}

			if (!hasAssociation) continue;

			for (ArmDefinition arm : arms) {
				if (arm == null) continue;

				String arm_oid = arm.getOid();
				if (arm_oid == null) throw new NullPointerException("An arm OID is NULL"); 

				List<ObservationSequence> oseqs = arm.getListOfObservationSequence();
				if (oseqs == null) continue;

				boolean armScopedAssociation = false;
				StandardAssignable start = null;
				for (ObservationSequence oseq : oseqs) {
					if (oseq == null) continue;
					start = oseq.getStart();
					if (armScopedAssociation) break;
					ObservationList ob_list = oseq.getObservationList();
					if (ob_list == null) continue;
					List<OidRef> ob_refs = ob_list.getListOfObservationRef();
					if (ob_refs == null) continue;
					for (OidRef ob_ref : ob_refs) {
						if (ob_ref == null) continue;
						String ob_ref_oid = ob_ref.getOidRef();
						if (window.equals(ob_ref_oid)) {
							armScopedAssociation = true;
							break;
						}
					}
				}

				if (armScopedAssociation) {
					List<StandardAssignable> tps = ob_timepoint_map.get(window);
					if (tps != null) {
						ElementaryDesign ed = new ElementaryDesign(oblock);
						
						ed.arm_oid = arm_oid;
						ed.output = referenced_output;
						ed.observation_oid = window;
						
						for (StandardAssignable tp : tps) ed.addTimepointExpression(tp);
						BinaryTree bt = tm.newInstance(start);
						ed.start_time_offset = Double.parseDouble(p.parse(ctx, bt).trim());
						
						protocol.addElementDesign(ed);
					}
				}
			}
		}

		if (!protocol.hasElementaryDesigns()) protocol = null;

		return protocol;
	}
	
	private Map<ObservationBlock, PharmMLRootType> createStandardMap() {
		Map<ObservationBlock, PharmMLRootType> standard_map = new HashMap<ObservationBlock, PharmMLRootType>(); 
		List<ObservationBlock> oblocks = lexer.getObservationBlocks();
		if (oblocks.isEmpty()) return standard_map;
		
		for (ObservationBlock block : oblocks) {
			if (block == null) continue;
		
			if (!isStructuredError(block.getObservationError())) 
				throw new IllegalStateException("Unable to build a protocol as error model is not a structured error model.");
			
			StructuredObsError soe = (StructuredObsError) block.getObservationError();
			standard_map.put(block, soe);
		}
		
		return standard_map;
	}
	
	private Map<String, List<String>> createWindowContainmentMap() {
		Map<String, List<String>> map = new HashMap<String, List<String>>();
		Observations observations = td.getObservations();
		if (observations == null) return map;
		
		if (observations.getListOfObservationsCombination() != null) {
			for (ObservationsCombination oc : observations.getListOfObservationsCombination()) {
				if (oc == null) continue;
				if (oc.getListOfObservations() == null) continue;
				String oid = oc.getOid(); 
				if (oid == null) throw new NullPointerException("An observation combination identifier is NULL");
				List<String> windows = new ArrayList<String>();
				for (SingleObservation so : oc.getListOfObservations()) {
					if (so == null) continue;
					if (so.getListOfObservationRef() == null) continue;
					for (OidRef ref : so.getListOfObservationRef()) {
						if (ref == null) continue;
						String inner_oid = ref.getOidRef();
						
						if (inner_oid == null) continue;
						if (!windows.contains(inner_oid)) windows.add(inner_oid);
					}
				}

				if (!windows.contains(oid)) windows.add(oid);
				map.put(oid, windows);
			}
		}
		
		if (observations.getListOfObservationsElements() != null) {
			for (PharmMLRootType o : observations.getListOfObservationsElements()) {
				if (isObservation(o)) {
					Observation ob = (Observation) o;
					String oid = ob.getOid();
					if (oid == null) throw new NullPointerException("An observation identifier is NULL");
					List<String> windows = new ArrayList<String>();
					windows.add(oid);
					map.put(oid, windows);
				}
			}
		}
		
		return map;
	}
	
	@Override
	public Administration getAdministration(String oid) {
		if (oid == null) return null;
		else return admin_map_.get(oid);
	}
	
	@Override
	public Map<String, Administration> getAdministrationMap() { return admin_map_; }
	
	@Override
	public List<Administration> getAdministrations() { return admins; }
	
	@Override
	public double getAdministrationStartTime(String admin_oid) {
		if (admin_oid == null) return 0.0;
		if (dose_time_raw.containsKey(admin_oid)) return dose_time_raw.get(admin_oid);
		else return 0.0;
	} 
	
	@Override
	public ArmDefinition getArm(String arm_oid) {
		if (arm_oid == null) return null;
		else if (arm_map.containsKey(arm_oid)) return arm_map.get(arm_oid);
		return null;
	}
	
	@Override
	public int getArmCount() { return arms.size(); }
	
	@Override
	public List<ArmDefinition> getArmMembership(Observation ob) {
		List<ArmDefinition> list = new ArrayList<ArmDefinition>();
		
		if (ob != null) {
			if (arm_membership_map.containsKey(ob)) {
				List<ArmDefinition> list2 = arm_membership_map.get(ob);
				if (list2 != null) list.addAll(list2);
			}
		}
		
		return list;
	}
	
	@Override
	public Map<Observation, List<ArmDefinition>> getArmMembershipMap() { return arm_membership_map; }
	
	@Override
	public List<ArmDefinition> getArms() { return arms; }
	
	@Override
	public int getArmSize(String oid) {
		if (oid == null) return 0;
		int size = 0;
		if (arm_size_map.containsKey(oid)) size = arm_size_map.get(oid);
		return size;
	}
	
	@Override
	public Map<SingleDesignSpace, Observation> getDesignSpaceObservationMap() { return designspace_2_window; }
	
	@Override
	public List<SingleDesignSpace> getDesignSpaces() { return design_spaces; }
	
	@Override
	public Map<ElementaryDesign, List<SingleDesignSpace>> getDesignSpaces(Protocol protocol) {
		if (protocol == null) throw new NullPointerException("Protocol is NULL");
		if (!protocol.hasElementaryDesigns()) throw new NullPointerException("Protocol has no elementary design");
		
		Map<String, List<String>> windows_map = createWindowContainmentMap();
		Map<ElementaryDesign, List<SingleDesignSpace>> map = new HashMap<ElementaryDesign, List<SingleDesignSpace>>();
		List<SingleDesignSpace> spaces = getDesignSpaces();
		
		for (ElementaryDesign ed : protocol.elementary_designs) {
			if (ed == null) continue;
			String window = ed.observation_oid;
			if (window == null) throw new NullPointerException("Elementary design Observation OID value was NULL");
			List<String> ds_referenced_windows = windows_map.get(window);
			List<SingleDesignSpace> ds = new ArrayList<SingleDesignSpace>();
			if (ds_referenced_windows == null) throw new NullPointerException("Unable create a list of desidn space referenced windows");
			for (SingleDesignSpace space : spaces) {
				if (space == null) continue;
				List<OidRef> windows = space.getListOfObservationRef();
				if (windows == null) continue;
				for (OidRef ref : windows) {
					if (ref == null) continue;
					String window_oid = ref.getOidRef();
					if (window_oid == null) continue;
					if (ds_referenced_windows.contains(window_oid)) {
						if (!ds.contains(space)) ds.add(space);
					}
				}
			}
			
			map.put(ed, ds);
		}
		
		return map;
	}

	@Override
	public SingleDesignSpace getDesignSpaceWithSamplingCountLimits() {
		for (SingleDesignSpace design_space : design_spaces) {
			if (design_space == null) continue;
			if (design_space.getNumberTimes() != null) return design_space;
		}
		
		return null;
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
	
	@Override
	public Map<String, PharmMLRootType> getDoseTargetMap() { return dose_target_map; }
	
	@Override
	public List<PharmMLRootType> getDoseTargets() {
		List<PharmMLRootType> targets = new ArrayList<PharmMLRootType>();
		for (PharmMLRootType target : dose_target_map.values()) {
			if (target == null) continue;
			if (!targets.contains(target)) targets.add(target);
		}
		return targets;
	}
	
	@Override
	public InterventionSequenceRef getInterventionSequenceRef(ArmDefinition arm) {
		if (arm == null) return null;
		String oid = arm.getOid();
		if (iseq_map.containsKey(oid)) return iseq_map.get(oid);
		else return null;
	}
	
	@Override
	public TrialDesign getModel() { return td; }
	
	@Override
	public String getName() { return "trial_design"; }
	
	@Override
	public Observation getObservation(ArmDefinition arm) {
		if (arm == null) return null;
		if (arm_2_observation_map.containsKey(arm)) {
			String ob_oid = arm_2_observation_map.get(arm);
			return this.getObservation(ob_oid);
		}
		
		return null;
	}

	@Override
	public Observation getObservation(OidRef ref) {
		if (ref == null) return null;
		else return getObservation(ref.getOidRef());
	}
	
	@Override
	public Observation getObservation(SingleDesignSpace space) {
		if (space == null) return null;
		else if (designspace_2_window.containsKey(space)) return designspace_2_window.get(space); 
		return null;
	}

	@Override
	public Observation getObservation(String oid) {
		if (oid == null) return null;
		else if (obs_map.containsKey(oid))  return obs_map.get(oid);
		else return null;
	}

	@Override
	public int getObservationIndex(OidRef ob_ref) {
		if (ob_ref == null) return -1;
		
		String oid = ob_ref.getOidRef();
		if (oid == null) return -1;
		
		int i = 0;
		for (Observation ob : obs) {
			if (oid.equals(ob.getOid())) return i;
			i++;
		}
		
		return -1;
	}
	
	@Override
	public List<Observation> getObservations() { return obs; }

	@Override
	public double getObservationStart(ArmDefinition arm) {
		double start = 0.0;
		if (arm != null) {
			if (arm_observation_start_map.containsKey(arm)) start = arm_observation_start_map.get(arm);
		}
		
		return start;
	}

	@Override
	public List<Protocol> getProtocols() { return protocols; }

	@Override
	public List<DerivativeVariable> getStateVariablesWithDosing() { throw new UnsupportedOperationException(); }

	@Override
	public List<String> getSymbolIds() { return new ArrayList<String>(); }

	@Override
	public boolean hasAdministration(String oid) {
		if (oid == null) return false;
		return admin_map_.containsKey(oid);
	}

	@Override
	public boolean hasDesignSpaceWithSamplingCountLimits() {
		for (SingleDesignSpace design_space : design_spaces) {
			if (design_space == null) continue;
			if (design_space.getNumberTimes() != null) return true;
		}
		
		return false;
	}

	@Override
	public boolean hasDosing() { return getStateVariablesWithDosing().isEmpty() == false; }

	@Override
	public boolean hasOccassions() { return td.getListOfOccasions().size() > 0; }

	@Override
	public boolean hasSymbolId(String name) { return false; }

	private void initArmMembershipMap() {
		if (obs.isEmpty()) return;
		
		for (Observation ob : obs) {
			if (ob == null) continue;
			if (!arm_membership_map.containsKey(ob)) arm_membership_map.put(ob, new ArrayList<ArmDefinition>());
		}
	}

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

	
	private void processDesignSpaces(List<SingleDesignSpace> spaces) {
		if (spaces == null) return;
		if (spaces.isEmpty()) return;
		
		for (SingleDesignSpace space : spaces) {
			if (space == null) continue;
			design_spaces.add(space);
			
			lexer.addStatement(space, tm.newInstance(space));
			lexer.updateNestedTrees();
			registerObservationWithDesignSpace(space);
		}
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
	
	private void registerObservationWithDesignSpace(SingleDesignSpace space) {
		if (space == null) return;
		List<OidRef> refs = space.getListOfObservationRef();
		if (refs == null) throw new IllegalStateException("A design space does not reference an observation window.");
		if (refs.isEmpty()) throw new IllegalStateException("A design space does not reference an observation window.");
		
		OidRef ref = refs.get(0);
		if (ref == null) throw new NullPointerException("A design space observation reference is NULL");
		Observation ob = getObservation(ref);
		if (ob == null) throw new IllegalStateException("An design space observation reference is not valid");
		
		designspace_2_window.put(space, ob);
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

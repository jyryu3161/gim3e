"""
GIM3E Example: Recon 2 Human Metabolic Model with Transcriptomics and Metabolomics Integration

This example demonstrates how to use GIM3E with the Recon 2 human metabolic model
to integrate transcriptomics (RNA-seq or microarray) and metabolomics data for
condition-specific metabolic modeling.

Author: Updated for modern COBRApy
Date: 2025
Requirements:
    - cobra >= 0.26.0
    - pandas
    - numpy
    - cplex or gurobi (recommended for MILP)
"""

import os
import pickle
import pandas as pd
import numpy as np
from copy import deepcopy

# Import GIM3E core functionality
from gim3e.core import gim3e

# Import COBRApy utilities
try:
    from cobra.io import read_sbml_model, load_model
    from cobra import Model, Reaction, Metabolite
except ImportError:
    print("Error: Please install modern COBRApy (>= 0.26.0)")
    print("  pip install cobra>=0.26.0")
    exit(1)


def load_recon2_model(model_path=None):
    """
    Load the Recon 2 human metabolic model.

    Parameters
    ----------
    model_path : str, optional
        Path to Recon 2 SBML file. If None, will try to download from BiGG Models.

    Returns
    -------
    cobra.Model
        Recon 2 metabolic model

    Notes
    -----
    Recon 2 can be obtained from:
    - BiGG Models: http://bigg.ucsd.edu/models/Recon2
    - VMH: https://www.vmh.life/

    For this example, we'll use Recon2.2 which is more curated.
    """
    if model_path and os.path.exists(model_path):
        print(f"Loading Recon 2 model from {model_path}...")
        model = read_sbml_model(model_path)
    else:
        # Try to load from BiGG Models database
        try:
            print("Attempting to load Recon2.2 from BiGG Models...")
            model = load_model("Recon2M2")  # Recon2M2 is a maintained version
            print(f"Successfully loaded {model.id}")
        except Exception as e:
            print(f"Could not load model: {e}")
            print("\nPlease download Recon 2 manually from:")
            print("  http://bigg.ucsd.edu/models/Recon2M2")
            print("Or provide a local SBML file path")

            # Create a minimal example model for demonstration
            print("\nCreating minimal example model for demonstration...")
            model = create_minimal_human_model()

    return model


def create_minimal_human_model():
    """
    Create a minimal human-like metabolic model for demonstration purposes.

    This is a toy model to demonstrate the GIM3E workflow when Recon 2 is not available.
    """
    model = Model("minimal_human_example")

    # Create compartments
    # Extracellular metabolites
    glc_e = Metabolite('glc__D_e', formula='C6H12O6', name='D-Glucose', compartment='e')
    o2_e = Metabolite('o2_e', formula='O2', name='Oxygen', compartment='e')
    co2_e = Metabolite('co2_e', formula='CO2', name='CO2', compartment='e')
    h2o_e = Metabolite('h2o_e', formula='H2O', name='Water', compartment='e')
    lac_e = Metabolite('lac__D_e', formula='C3H6O3', name='D-Lactate', compartment='e')
    pyr_e = Metabolite('pyr_e', formula='C3H4O3', name='Pyruvate', compartment='e')
    gln_e = Metabolite('gln__L_e', formula='C5H10N2O3', name='L-Glutamine', compartment='e')
    glu_e = Metabolite('glu__L_e', formula='C5H9NO4', name='L-Glutamate', compartment='e')

    # Cytoplasmic metabolites
    glc_c = Metabolite('glc__D_c', formula='C6H12O6', name='D-Glucose', compartment='c')
    g6p_c = Metabolite('g6p_c', formula='C6H13O9P', name='Glucose-6-phosphate', compartment='c')
    f6p_c = Metabolite('f6p_c', formula='C6H13O9P', name='Fructose-6-phosphate', compartment='c')
    pyr_c = Metabolite('pyr_c', formula='C3H4O3', name='Pyruvate', compartment='c')
    lac_c = Metabolite('lac__D_c', formula='C3H6O3', name='D-Lactate', compartment='c')
    atp_c = Metabolite('atp_c', formula='C10H16N5O13P3', name='ATP', compartment='c')
    adp_c = Metabolite('adp_c', formula='C10H15N5O10P2', name='ADP', compartment='c')
    nad_c = Metabolite('nad_c', formula='C21H27N7O14P2', name='NAD', compartment='c')
    nadh_c = Metabolite('nadh_c', formula='C21H28N7O14P2', name='NADH', compartment='c')

    # Mitochondrial metabolites
    pyr_m = Metabolite('pyr_m', formula='C3H4O3', name='Pyruvate', compartment='m')
    accoa_m = Metabolite('accoa_m', formula='C23H38N7O17P3S', name='Acetyl-CoA', compartment='m')

    # Exchange reactions
    ex_glc = Reaction('EX_glc__D_e')
    ex_glc.name = 'Glucose exchange'
    ex_glc.add_metabolites({glc_e: -1})
    ex_glc.lower_bound = -10  # uptake
    ex_glc.upper_bound = 0

    ex_o2 = Reaction('EX_o2_e')
    ex_o2.name = 'Oxygen exchange'
    ex_o2.add_metabolites({o2_e: -1})
    ex_o2.lower_bound = -20
    ex_o2.upper_bound = 0

    ex_lac = Reaction('EX_lac__D_e')
    ex_lac.name = 'Lactate exchange'
    ex_lac.add_metabolites({lac_e: -1})
    ex_lac.lower_bound = 0
    ex_lac.upper_bound = 1000

    # Transport reactions
    glc_t = Reaction('GLCt1')
    glc_t.name = 'Glucose transport'
    glc_t.add_metabolites({glc_e: -1, glc_c: 1})
    glc_t.gene_reaction_rule = 'SLC2A1'  # GLUT1

    # Glycolysis
    hk = Reaction('HEX1')
    hk.name = 'Hexokinase'
    hk.add_metabolites({glc_c: -1, atp_c: -1, g6p_c: 1, adp_c: 1})
    hk.gene_reaction_rule = 'HK1 or HK2'

    pgi = Reaction('PGI')
    pgi.name = 'Phosphoglucose isomerase'
    pgi.add_metabolites({g6p_c: -1, f6p_c: 1})
    pgi.gene_reaction_rule = 'GPI'

    # Simplified glycolysis to pyruvate
    glyc = Reaction('GLYC_simplified')
    glyc.name = 'Simplified glycolysis'
    glyc.add_metabolites({f6p_c: -1, atp_c: 2, nad_c: -1, nadh_c: 1, pyr_c: 2})
    glyc.gene_reaction_rule = 'PFKM and PKM'

    # Lactate dehydrogenase
    ldh = Reaction('LDH_D')
    ldh.name = 'Lactate dehydrogenase'
    ldh.add_metabolites({pyr_c: -1, nadh_c: -1, lac_c: 1, nad_c: 1})
    ldh.gene_reaction_rule = 'LDHA or LDHB'

    lac_t = Reaction('LACt')
    lac_t.name = 'Lactate transport'
    lac_t.add_metabolites({lac_c: -1, lac_e: 1})
    lac_t.gene_reaction_rule = 'SLC16A1'  # MCT1

    # ATP maintenance
    atpm = Reaction('ATPM')
    atpm.name = 'ATP maintenance'
    atpm.add_metabolites({atp_c: -1, adp_c: 1})
    atpm.lower_bound = 1.0

    # Biomass reaction (simplified)
    biomass = Reaction('BIOMASS')
    biomass.name = 'Biomass production'
    biomass.add_metabolites({
        atp_c: -50,
        adp_c: 50,
        g6p_c: -5,
    })
    biomass.lower_bound = 0
    biomass.upper_bound = 1000

    # Add all reactions to model
    model.add_reactions([
        ex_glc, ex_o2, ex_lac,
        glc_t, hk, pgi, glyc, ldh, lac_t,
        atpm, biomass
    ])

    # Set objective
    model.objective = 'BIOMASS'

    print(f"Created minimal model with {len(model.reactions)} reactions and {len(model.metabolites)} metabolites")

    return model


def generate_example_transcriptomics_data(model, condition="cancer"):
    """
    Generate example transcriptomics data for demonstration.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    condition : str
        Condition to simulate ('cancer', 'normal', 'hypoxia')

    Returns
    -------
    dict
        Gene expression dictionary {gene_id: expression_value}

    Notes
    -----
    In real applications, this data would come from:
    - RNA-seq (TPM, FPKM, or counts normalized)
    - Microarray (log2 intensity values)
    - Proteomics (protein abundance)

    Expression values can be:
    - Log2 fold change vs control
    - Log2 intensity/abundance
    - -log10(p-value) for presence/absence
    """
    expression_dict = {}

    # Get all genes from the model
    all_genes = [gene.id for gene in model.genes]

    if condition == "cancer":
        # Cancer cells often have high glycolysis (Warburg effect)
        # and glutaminolysis
        expression_patterns = {
            # High glycolysis genes
            'HK2': 8.5,      # Hexokinase 2 (high in cancer)
            'PFKM': 7.8,     # Phosphofructokinase
            'PKM': 8.2,      # Pyruvate kinase M2
            'LDHA': 9.0,     # Lactate dehydrogenase A (Warburg)
            'SLC2A1': 8.5,   # GLUT1 (high glucose uptake)
            'SLC16A1': 7.5,  # MCT1 (lactate export)

            # Moderate/low expression
            'HK1': 5.0,      # Hexokinase 1 (lower in cancer)
            'LDHB': 4.5,     # LDH B (lower in cancer)
            'GPI': 7.0,      # Phosphoglucose isomerase
        }
    elif condition == "normal":
        # Normal cells have balanced metabolism
        expression_patterns = {
            'HK1': 7.5,
            'HK2': 5.0,
            'PFKM': 6.5,
            'PKM': 6.5,
            'LDHA': 6.0,
            'LDHB': 7.0,
            'SLC2A1': 6.5,
            'SLC16A1': 6.0,
            'GPI': 7.0,
        }
    elif condition == "hypoxia":
        # Hypoxic cells upregulate glycolysis and HIF targets
        expression_patterns = {
            'HK2': 9.5,      # HIF target
            'PFKM': 8.5,
            'PKM': 8.8,
            'LDHA': 9.5,     # HIF target
            'SLC2A1': 9.0,   # HIF target
            'SLC16A1': 8.0,
            'HK1': 4.5,
            'LDHB': 4.0,
            'GPI': 7.5,
        }

    # Add baseline expression for all genes
    for gene in all_genes:
        if gene in expression_patterns:
            expression_dict[gene] = expression_patterns[gene]
        else:
            # Random baseline expression for other genes
            expression_dict[gene] = np.random.normal(6.0, 1.0)

    return expression_dict


def generate_example_metabolomics_data(model, condition="cancer"):
    """
    Generate example metabolomics data for demonstration.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    condition : str
        Condition to simulate

    Returns
    -------
    list
        List of detected metabolite IDs

    Notes
    -----
    In real applications, metabolomics data comes from:
    - LC-MS/MS (liquid chromatography-mass spectrometry)
    - GC-MS (gas chromatography-mass spectrometry)
    - NMR (nuclear magnetic resonance)

    The metabolite list represents metabolites detected above
    background levels, indicating they are actively used in metabolism.
    """
    # Get cytoplasmic metabolites (most common in metabolomics)
    cytoplasmic_mets = [m.id for m in model.metabolites if m.compartment == 'c']

    if condition == "cancer":
        # Cancer cells accumulate lactate, consume glucose and glutamine
        detected_metabolites = [
            'glc__D_c',      # Glucose (consumed)
            'g6p_c',         # Glucose-6-phosphate
            'f6p_c',         # Fructose-6-phosphate
            'pyr_c',         # Pyruvate
            'lac__D_c',      # Lactate (accumulated)
            'atp_c',         # ATP
            'adp_c',         # ADP
            'nad_c',         # NAD+
            'nadh_c',        # NADH
        ]
    else:
        # Normal metabolite detection
        detected_metabolites = [
            'glc__D_c',
            'g6p_c',
            'pyr_c',
            'atp_c',
            'adp_c',
        ]

    # Filter to only metabolites that exist in the model
    detected_metabolites = [m for m in detected_metabolites
                          if m in [met.id for met in model.metabolites]]

    return detected_metabolites


def setup_culture_medium(model, medium_type="RPMI"):
    """
    Set up culture medium conditions.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    medium_type : str
        Type of culture medium ('RPMI', 'DMEM', 'MEM')

    Notes
    -----
    Common culture media for human cells:
    - RPMI-1640: 11 mM glucose, 2 mM glutamine
    - DMEM: 25 mM glucose, 4 mM glutamine (high glucose variant)
    - MEM: 5.5 mM glucose, 2 mM glutamine
    """
    # Set all exchanges to 0 first
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 1000

    # Configure based on medium type
    if medium_type == "RPMI":
        glucose_uptake = 11.0  # mM
        glutamine_uptake = 2.0
    elif medium_type == "DMEM":
        glucose_uptake = 25.0
        glutamine_uptake = 4.0
    else:  # MEM or default
        glucose_uptake = 5.5
        glutamine_uptake = 2.0

    # Set nutrient availability
    try:
        model.reactions.get_by_id('EX_glc__D_e').lower_bound = -glucose_uptake
        model.reactions.get_by_id('EX_o2_e').lower_bound = -20.0
    except KeyError:
        print("Warning: Some exchange reactions not found in model")

    print(f"Set up {medium_type} medium: {glucose_uptake} mM glucose")


def run_gim3e_analysis(model, expression_dict, metabolite_list,
                       condition_name="cancer", solver='glpk'):
    """
    Run GIM3E analysis with transcriptomics and metabolomics integration.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model
    expression_dict : dict
        Gene expression data
    metabolite_list : list
        Detected metabolites
    condition_name : str
        Name of the condition
    solver : str
        Solver to use ('cplex', 'gurobi', 'glpk')

    Returns
    -------
    tuple
        (gim3e_model, FVA_results, total_penalty)
    """
    print(f"\n{'='*70}")
    print(f"Running GIM3E Analysis for {condition_name} condition")
    print(f"{'='*70}")

    # Calculate expression threshold (e.g., median expression)
    expression_values = list(expression_dict.values())
    expression_threshold = np.median(expression_values)

    print(f"\nExpression statistics:")
    print(f"  Total genes: {len(expression_dict)}")
    print(f"  Mean expression: {np.mean(expression_values):.2f}")
    print(f"  Median expression: {expression_threshold:.2f}")
    print(f"  Genes below threshold: {sum(1 for v in expression_values if v < expression_threshold)}")

    print(f"\nMetabolomics data:")
    print(f"  Detected metabolites: {len(metabolite_list)}")
    print(f"  Metabolites: {', '.join(metabolite_list[:5])}...")

    # GIM3E parameters
    fraction_growth = 0.9         # Require 90% of optimal growth
    relative_penalty_bound = 1.05  # Allow 5% penalty increase
    solver_tolerance = 1e-7

    print(f"\nGIM3E parameters:")
    print(f"  Expression threshold: {expression_threshold:.2f}")
    print(f"  Minimum growth: {fraction_growth*100:.0f}% of optimal")
    print(f"  Penalty bound: {relative_penalty_bound}")
    print(f"  MILP formulation: True")
    print(f"  Solver: {solver}")

    try:
        # Run GIM3E
        print(f"\nRunning GIM3E optimization...")
        gim3e_model, FVA_results, total_penalty = gim3e.gim3e(
            model,
            expression_dict=expression_dict,
            expression_threshold=expression_threshold,
            metabolite_list=metabolite_list,
            fraction_growth=fraction_growth,
            relative_penalty_bound=relative_penalty_bound,
            solver_tolerance=solver_tolerance,
            metabolite_flux_requirement=True,
            monitor_all_cellular_metabolites=False,
            MILP_formulation=True,
            run_FVA=True,
            reduce_model=False,
            trim_model=False,
            solver=solver
        )

        print(f"\n{'='*70}")
        print(f"GIM3E Analysis Complete!")
        print(f"{'='*70}")
        print(f"Total penalty score: {total_penalty:.4f}")
        print(f"FVA performed on {len(FVA_results)} reactions")

        return gim3e_model, FVA_results, total_penalty

    except Exception as e:
        print(f"\nError running GIM3E: {e}")
        print("\nTroubleshooting tips:")
        print("1. Ensure you have a compatible solver installed (cplex, gurobi, or glpk)")
        print("2. For MILP problems, CPLEX or Gurobi strongly recommended")
        print("3. Check that model is feasible: model.optimize()")
        print("4. Verify metabolite IDs match model metabolites")
        raise


def analyze_results(model, FVA_results, condition_name):
    """
    Analyze and visualize GIM3E results.

    Parameters
    ----------
    model : cobra.Model
        GIM3E processed model
    FVA_results : dict
        Flux variability analysis results
    condition_name : str
        Condition name for labeling
    """
    print(f"\n{'='*70}")
    print(f"Analysis Results - {condition_name}")
    print(f"{'='*70}")

    # Analyze active reactions
    active_reactions = []
    for rxn_id, fva in FVA_results.items():
        if fva['maximum'] is not None and abs(fva['maximum']) > 1e-6:
            active_reactions.append(rxn_id)

    print(f"\nActive reactions: {len(active_reactions)} / {len(FVA_results)}")

    # Find key metabolic fluxes
    print(f"\nKey metabolic fluxes:")
    key_reactions = ['HEX1', 'PGI', 'LDH_D', 'GLYC_simplified', 'BIOMASS', 'ATPM']

    for rxn_id in key_reactions:
        if rxn_id in FVA_results:
            fva = FVA_results[rxn_id]
            rxn = model.reactions.get_by_id(rxn_id)
            print(f"  {rxn.name} ({rxn_id}):")
            print(f"    Min: {fva['minimum']:.4f}, Max: {fva['maximum']:.4f}")

    # Analyze turnover metabolites
    tm_reactions = [r for r in model.reactions if r.id.startswith('TMS_')]
    if tm_reactions:
        print(f"\nTurnover metabolites (detected in metabolomics):")
        for rxn in tm_reactions[:10]:  # Show first 10
            met_id = rxn.id.replace('TMS_', '')
            if met_id in FVA_results:
                fva = FVA_results[met_id]
                print(f"  {met_id}: flux range [{fva['minimum']:.4f}, {fva['maximum']:.4f}]")


def save_results(model, FVA_results, expression_dict, metabolite_list,
                condition_name, output_dir="output"):
    """
    Save GIM3E results to files.

    Parameters
    ----------
    model : cobra.Model
        GIM3E model
    FVA_results : dict
        FVA results
    expression_dict : dict
        Expression data
    metabolite_list : list
        Metabolomics data
    condition_name : str
        Condition name
    output_dir : str
        Output directory
    """
    os.makedirs(output_dir, exist_ok=True)

    # Save model
    model_file = os.path.join(output_dir, f"gim3e_model_{condition_name}.pickle")
    with open(model_file, 'wb') as f:
        pickle.dump(model, f)
    print(f"\nSaved model to: {model_file}")

    # Save FVA results
    fva_file = os.path.join(output_dir, f"FVA_results_{condition_name}.pickle")
    with open(fva_file, 'wb') as f:
        pickle.dump(FVA_results, f)
    print(f"Saved FVA results to: {fva_file}")

    # Save to CSV for easy viewing
    fva_df = pd.DataFrame.from_dict(FVA_results, orient='index')
    csv_file = os.path.join(output_dir, f"FVA_results_{condition_name}.csv")
    fva_df.to_csv(csv_file)
    print(f"Saved FVA results (CSV) to: {csv_file}")


def main():
    """
    Main function to demonstrate GIM3E with Recon 2 and multi-omics integration.
    """
    print("="*70)
    print("GIM3E: Recon 2 Transcriptomics-Metabolomics Integration Example")
    print("="*70)

    # Configuration
    CONDITION = "cancer"  # Options: 'cancer', 'normal', 'hypoxia'
    SOLVER = 'glpk'       # Options: 'cplex', 'gurobi', 'glpk'
    MEDIUM = "RPMI"       # Options: 'RPMI', 'DMEM', 'MEM'

    print(f"\nConfiguration:")
    print(f"  Condition: {CONDITION}")
    print(f"  Solver: {SOLVER}")
    print(f"  Culture medium: {MEDIUM}")

    # Step 1: Load Recon 2 model
    print(f"\n{'='*70}")
    print("Step 1: Loading Recon 2 Model")
    print(f"{'='*70}")
    model = load_recon2_model()

    print(f"\nModel statistics:")
    print(f"  Reactions: {len(model.reactions)}")
    print(f"  Metabolites: {len(model.metabolites)}")
    print(f"  Genes: {len(model.genes)}")

    # Step 2: Set up culture medium
    print(f"\n{'='*70}")
    print("Step 2: Setting Up Culture Medium")
    print(f"{'='*70}")
    setup_culture_medium(model, MEDIUM)

    # Test model optimization
    print("\nTesting model optimization...")
    solution = model.optimize()
    print(f"Model optimization status: {solution.status}")
    if solution.status == 'optimal':
        print(f"Objective value (growth): {solution.objective_value:.4f}")

    # Step 3: Generate/load omics data
    print(f"\n{'='*70}")
    print("Step 3: Preparing Omics Data")
    print(f"{'='*70}")

    # Generate example transcriptomics data
    expression_dict = generate_example_transcriptomics_data(model, CONDITION)
    print(f"Generated transcriptomics data for {len(expression_dict)} genes")

    # Generate example metabolomics data
    metabolite_list = generate_example_metabolomics_data(model, CONDITION)
    print(f"Generated metabolomics data: {len(metabolite_list)} metabolites detected")

    # Step 4: Run GIM3E
    try:
        gim3e_model, FVA_results, total_penalty = run_gim3e_analysis(
            model,
            expression_dict,
            metabolite_list,
            CONDITION,
            solver=SOLVER
        )

        # Step 5: Analyze results
        analyze_results(gim3e_model, FVA_results, CONDITION)

        # Step 6: Save results
        save_results(gim3e_model, FVA_results, expression_dict,
                    metabolite_list, CONDITION)

        print("\n" + "="*70)
        print("Analysis Complete!")
        print("="*70)
        print("\nNext steps:")
        print("1. Examine FVA results to identify condition-specific active pathways")
        print("2. Compare multiple conditions (e.g., cancer vs normal)")
        print("3. Perform metabolite knockout analysis to find essential metabolites")
        print("4. Validate predictions with experimental data")

    except Exception as e:
        print(f"\nAnalysis failed: {e}")
        print("\nNote: This example requires a working solver.")
        print("For production use, install CPLEX or Gurobi.")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

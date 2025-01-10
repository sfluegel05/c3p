"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid conjugated to a functional group that adds hydrophilicity or charge.
    Conjugating groups include glycine, taurine, sulfuric acid, glucuronic acid, glucose, amino acids, and coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdqueries
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone (four fused rings: three 6-membered and one 5-membered)
    ri = mol.GetRingInfo()
    if not ri.IsInitialized():
        return False, "No ring information available"

    num_rings = ri.NumRings()
    if num_rings < 4:
        return False, "Less than 4 rings found, not a steroid backbone"

    # Get ring sizes
    ring_sizes = [len(ri.AtomRings()[i]) for i in range(num_rings)]
    # Check for three 6-membered rings and one 5-membered ring
    if ring_sizes.count(6) < 3 or ring_sizes.count(5) < 1:
        return False, "Does not have the required steroid ring sizes (three 6-membered and one 5-membered rings)"
    
    # Check for fused rings
    fused = True
    atom_rings = ri.AtomRings()
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            if set(atom_rings[i]) & set(atom_rings[j]):
                continue
            else:
                fused = False
                break
    if not fused:
        return False, "Rings are not fused as in steroid backbone"
    
    # Look for carboxylic acid or conjugated side chain
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O-]")  # Deprotonated
    carboxylic_acid_protonated = Chem.MolFromSmarts("C(=O)O")  # Protonated
    if mol.HasSubstructMatch(carboxylic_acid) or mol.HasSubstructMatch(carboxylic_acid_protonated):
        # Check for conjugation with specified groups
        # Define patterns for conjugating groups
        glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
        taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)O")
        sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")
        glucuronic_acid_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
        glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
        coenzymeA_pattern = Chem.MolFromSmarts("NC(=O)CNC(=O)CC(N)C(=O)NCCS")
        amino_acid_pattern = Chem.MolFromSmarts("NCC(=O)O")

        # Check for amide bond to amino acids (including glycine and taurine)
        amide_bond = Chem.MolFromSmarts("C(=O)N")
        ester_bond = Chem.MolFromSmarts("C(=O)O")
        sulfonate_bond = Chem.MolFromSmarts("S(=O)(=O)O")
        
        conjugated = False
        reason = ""
        # Check for glycine conjugation
        if mol.HasSubstructMatch(glycine_pattern):
            conjugated = True
            reason = "Conjugated with glycine"
        # Check for taurine conjugation
        elif mol.HasSubstructMatch(taurine_pattern):
            conjugated = True
            reason = "Conjugated with taurine"
        # Check for sulfate group
        elif mol.HasSubstructMatch(sulfate_pattern):
            conjugated = True
            reason = "Contains sulfate group"
        # Check for glucuronic acid conjugation
        elif mol.HasSubstructMatch(glucuronic_acid_pattern):
            conjugated = True
            reason = "Conjugated with glucuronic acid"
        # Check for glucose conjugation
        elif mol.HasSubstructMatch(glucose_pattern):
            conjugated = True
            reason = "Conjugated with glucose"
        # Check for coenzyme A conjugation
        elif mol.HasSubstructMatch(coenzymeA_pattern):
            conjugated = True
            reason = "Conjugated with coenzyme A"
        # Check for other amino acid conjugations via amide bond
        elif mol.HasSubstructMatch(amide_bond) and mol.HasSubstructMatch(amino_acid_pattern):
            conjugated = True
            reason = "Conjugated with an amino acid"
        else:
            conjugated = False
            reason = "No conjugating groups found"
        
        if conjugated:
            return True, reason
        else:
            return False, reason
    else:
        return False, "No carboxylic acid or conjugated side chain found"

__metadata__ = {   'chemical_class': {   'name': 'bile acid conjugate',
                              'definition': 'Any bile acid conjugated to a '
                                            'functional group that gives additional '
                                            'hydrophilicity or charge to the molecule. '
                                            'Molecules used for conjugation are: glycine, '
                                            'taurine (and other amino acids); sulfuric '
                                            'acid (for which the term ''sulfate'' may be '
                                            'used); glucuronic acid (for which the term '
                                            ''''glucuronate'' may be used); glucose and '
                                            'other uncharged sugars; and coenzyme A.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None}
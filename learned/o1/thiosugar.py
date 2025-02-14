"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for pyranose ring with O or S in ring position
    pyranose_pattern = Chem.MolFromSmarts("C1[C@H]([OH])[C@@H](O)[C@H](O)[C@@H](O)[O,S]1")
    furanose_pattern = Chem.MolFromSmarts("C1[C@H]([OH])[C@@H](O)[C@@H](O)[O,S]1")
    
    # Check if molecule contains a sugar ring with possible sulfur substitution
    sugar_rings = mol.GetSubstructMatches(pyranose_pattern) + mol.GetSubstructMatches(furanose_pattern)
    if not sugar_rings:
        return False, "No sugar ring (pyranose or furanose) found"
    
    # Examine each sugar ring
    for ring in sugar_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Check if heteroatom in ring is sulfur
        heteroatom = mol.GetAtomWithIdx(ring[-1])  # Last atom in pattern is [O,S]1
        if heteroatom.GetAtomicNum() == 16:
            return True, "Contains sugar ring with sulfur in ring (thiosugar)"
        
        # Check for sulfur substitutions at hydroxyl groups
        sulfur_substituted = False
        for atom in ring_atoms[:-1]:  # Exclude heteroatom
            if atom.GetSymbol() != 'C':
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' or neighbor.GetSymbol() == 'S':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.Single:
                        if neighbor.GetAtomicNum() == 16:
                            sulfur_substituted = True
                            break
            if sulfur_substituted:
                break
        if sulfur_substituted:
            return True, "Contains sugar ring with sulfur substitution at hydroxyl group (thiosugar)"
    
    return False, "No sulfur substitutions found in sugar ring"

__metadata__ = {   'chemical_class': {   'name': 'thiosugar',
                                         'definition': 'A carbohydrate derivative in which one or more of the oxygens or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.'},
                   'config': {   'llm_model_name': 'lbl/claude-sonnet',
                                 'f1_threshold': 0.8,
                                 'max_attempts': 5,
                                 'max_positive_instances': None,
                                 'max_positive_to_test': None,
                                 'max_negative_to_test': None,
                                 'max_positive_in_prompt': 50,
                                 'max_negative_in_prompt': 20,
                                 'max_instances_in_prompt': 100,
                                 'test_proportion': 0.1}}
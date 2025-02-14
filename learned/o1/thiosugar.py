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

    # Check if molecule contains sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found in molecule"

    # Define SMARTS patterns for sugar rings (allowing for sulfur in ring)
    pyranose_pattern = Chem.MolFromSmarts("[C@H]1([O,S])[C@@H]([O,S])[C@H]([O,S])[C@@H]([O,S])[C@H]1")
    furanose_pattern = Chem.MolFromSmarts("[C@H]1([O,S])[C@@H]([O,S])[C@H]([O,S])[C@H]1")

    # Check if molecule contains sugar rings
    sugar_ring_found = False
    thio_substitution_found = False

    ring_matches = mol.GetSubstructMatches(pyranose_pattern) + mol.GetSubstructMatches(furanose_pattern)
    if ring_matches:
        sugar_ring_found = True
        # For each ring, check for sulfur substitutions
        for match in ring_matches:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]

            # Check for sulfur in ring positions
            for atom in ring_atoms:
                if atom.GetAtomicNum() == 16:  # Sulfur atom
                    thio_substitution_found = True
                    break  # No need to check further

            if thio_substitution_found:
                break  # Found at least one thiosugar ring

            # Check for sulfur substituents attached to ring carbons
            for atom in ring_atoms:
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 16:  # Sulfur atom
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                # Sulfur attached to ring carbon, replacing hydroxyl group
                                thio_substitution_found = True
                                break
                    if thio_substitution_found:
                        break
            if thio_substitution_found:
                break

    if not sugar_ring_found:
        return False, "No sugar ring (pyranose or furanose) found in molecule"

    if thio_substitution_found:
        return True, "Contains sugar ring with sulfur substitution"

    return False, "No sulfur substitutions found in sugar ring"

__metadata__ = {
    'chemical_class': {
        'name': 'thiosugar',
        'definition': 'A carbohydrate derivative in which one or more of the oxygens or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.'
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}
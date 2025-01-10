"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

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

    # Define a general SMARTS pattern for monosaccharide units (pyranose and furanose rings)
    sugar_pattern = Chem.MolFromSmarts("""
    [
        # Six-membered ring sugars (pyranose)
        [$([C;H1,H2,H3][O;H0])]1[C,O;R][C,O;R][C,O;R][C,O;R][C,O;R]1,
        # Five-membered ring sugars (furanose)
        [$([C;H1,H2,H3][O;H0])]1[C,O;R][C,O;R][C,O;R][C,O;R]1
    ]
    """)

    # Find sugar units in the molecule
    matches = mol.GetSubstructMatches(sugar_pattern)
    if not matches:
        return False, "No sugar units found in molecule"

    # For each sugar unit, check if any oxygen or hydroxyl group is replaced by sulfur or -SR
    thio_substitution_found = False
    for match in matches:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        for atom in ring_atoms:
            atom_num = atom.GetAtomicNum()
            # Check for sulfur atoms in place of oxygen in ring
            if atom_num == 16:
                thio_substitution_found = True
                break
            # Check for sulfur substituents replacing hydroxyl groups
            elif atom_num == 6:  # Carbon atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in match:
                        continue  # Skip atoms in the ring
                    if neighbor.GetAtomicNum() == 16:  # Sulfur atom
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            # Sulfur attached to ring carbon, possibly replacing hydroxyl
                            thio_substitution_found = True
                            break
                if thio_substitution_found:
                    break
        if thio_substitution_found:
            break

    if thio_substitution_found:
        return True, "Contains sugar unit with sulfur substitution of oxygen or hydroxyl group"

    return False, "No sulfur substitutions found in sugar units"

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
"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose must have the specific stereochemistry at the C5 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly six carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"
    
    # Hexoses typically have 5 OH groups (one on each non-anomeric carbon)
    oh_groups = mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]'))
    if len(oh_groups) < 4:
        return False, f"Expected at least 4 OH groups, found {len(oh_groups)}"
    
    # Define enhanced patterns to check for furanose and pyranose including stereochemistry
    furanose_pattern = Chem.MolFromSmarts('C1[C@H](O)[C@H](O)[C@@H](O)C1O')
    pyranose_pattern = Chem.MolFromSmarts('C1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C1O')
    
    # Specifically look for the D stereochemistry at the C5 position
    d_chiral_pattern = Chem.MolFromSmarts('[C@H](O)[C@H](O)[C@@H](CO)')

    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        # Find D-configuration matches
        match_positions = mol.GetSubstructMatches(d_chiral_pattern)
        if match_positions:
            for pos in match_positions:
                atom_idx = pos[0]
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    return True, "Identified as D-hexose with correct stereochemistry at C5"
    
    return False, "No matching D-hexose ring pattern with correct 5th carbon stereochemistry found"
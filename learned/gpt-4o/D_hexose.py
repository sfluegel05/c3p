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
    A D-hexose must have specific stereochemistry at the C5 position.

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

    # Hexoses typically have exactly 6 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"
    
    # Hexoses typically have several OH groups
    oh_groups = mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]'))
    if len(oh_groups) < 4:
        return False, f"Expected at least 4 OH groups, found {len(oh_groups)}"

    # D-hexose specific patterns, checking furanose and pyranose configurations
    pyranose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)C1')
    furanose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)C1')

    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        # Check the C5 position specifically has D-configuration
        c5_d_chiral_pattern = Chem.MolFromSmarts('[C@H](O)C[CH]1')
        
        # Ensure there are accompanying chiral centers as well
        if mol.HasSubstructMatch(c5_d_chiral_pattern):
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            if len(chiral_centers) >= 4:
                return True, "Identified as D-hexose with correct stereochemistry at C5"

    return False, "No matching D-hexose ring pattern with correct 5th carbon stereochemistry found"
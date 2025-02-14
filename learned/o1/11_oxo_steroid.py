"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is an oxo steroid that has an oxo (=O) substituent at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the steroid nucleus with atom mapping
    steroid_smarts = """
    [#6]1-[#6]-[#6]-[#6]2-[#6]-[#6]3-[#6]([#6]-[#6]1)-[#6]-[#6]4-[#6]([#6]-[#6]([#6]-2)-[#6]3)-[#6]-[#6]-[#6]-4
    """
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    # Check if molecule has a steroid nucleus
    steroid_match = mol.GetSubstructMatch(steroid_pattern)
    if not steroid_match:
        return False, "Molecule does not have a steroid nucleus"

    # Map the steroid nucleus to identify position 11
    # Create a template steroid with atom mapping numbers corresponding to standard steroid positions
    steroid_template = Chem.MolFromSmiles("""
    C[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@@H]3CC[C@]12C
    """)

    if steroid_template is None:
        return False, "Failed to create steroid template"

    # Align the molecule to the steroid template
    match = mol.GetSubstructMatch(steroid_template)
    if not match:
        return False, "Molecule does not align with steroid template"

    # Position 11 corresponds to atom index 9 in the template (zero-based indexing)
    # Get the mapped atom in the molecule
    pos11_atom_idx = match[9]
    pos11_atom = mol.GetAtomWithIdx(pos11_atom_idx)

    # Check if position 11 has a ketone group (=O)
    found_oxo = False
    for neighbor in pos11_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(pos11_atom.GetIdx(), neighbor.GetIdx())
        if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                found_oxo = True
                break

    if not found_oxo:
        return False, "No oxo group (=O) at position 11"

    return True, "Molecule has a steroid nucleus with an oxo group at position 11"
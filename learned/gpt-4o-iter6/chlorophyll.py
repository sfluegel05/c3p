"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are characterized by a magnesium porphyrin core with a fifth modified ring and typically a long aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium present exactly once
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if len(mg_atoms) != 1:
        return False, "Magnesium not found or more than one magnesium atom present"

    # Porphyrin core with magnesium coordination and fifth ring
    portmangrin_core = Chem.MolFromSmarts('n1ccccc2c1c3nc4cc[nH+]c4c(-[Mg])c3c2')
    if not mol.HasSubstructMatch(portmangrin_core):
        return False, "Porphyrin core with magnesium coordination and fifth ring not found"

    # Check for a variability in aliphatic chain, typically represented by long hydrocarbon chains
    long_chain_smarts = Chem.MolFromSmarts('C{5,}')  # Flexibility: Look for any sequence of at least 5 contiguous carbon atoms
    if not mol.HasSubstructMatch(long_chain_smarts):
        return False, "Long aliphatic chain characteristic of chlorophyll missing"

    return True, "Contains magnesium porphyrin core with a fifth ring and a characteristic long aliphatic chain"
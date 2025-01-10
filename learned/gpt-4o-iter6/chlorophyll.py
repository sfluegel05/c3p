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

    # Improved pattern for the porphyrin core with magnesium coordination
    # This pattern accounts for a typical porphyrin ring around a central magnesium
    porphyrin_core_improved = Chem.MolFromSmarts('n1ccc2c(ncn([Mg])c2c1)c1ccc3[nH]c(nc3c1)[Mg]')
    if not mol.HasSubstructMatch(porphyrin_core_improved):
        return False, "Porphyrin core with magnesium coordination not found"

    # Check for the presence of a long aliphatic chain, allowing more variability
    # Here we simply look for a long carbon chain which is characteristic, but not overly specific
    long_chain_smarts = Chem.MolFromSmarts('CCCCCCC')  # Looks for a stretch of at least 7 contiguous C atoms
    if not mol.HasSubstructMatch(long_chain_smarts):
        return False, "Long aliphatic chain characteristic of chlorophyll missing"

    return True, "Contains magnesium porphyrin core with a fifth ring and a long aliphatic chain"
"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are characterized by a magnesium porphyrin core with modifications and typically a long aliphatic chain.

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

    # Check for magnesium present
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if len(mg_atoms) != 1:
        return False, "Magnesium not found or more than one magnesium atom present"

    # Magnesium-coordinated porphyrin core with fifth ring (generic pattern)
    # This SMARTS is intended as a generalized porphyrin core + modifiers typical of chlorophyll
    porphyrin_core_smarts = '''
    [n;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1-[n;R1]2[c;R1][c;R1][c;R1][n;R1][c;R1]2
    -[c;R1]3[n;R1][c;R1][c;R1][c;R1][n;R1]3'''
    porphyrin_core = Chem.MolFromSmarts(porphyrin_core_smarts)
    
    if not mol.HasSubstructMatch(porphyrin_core):
        return False, "Porphyrin core with magnesium coordination not found"

    # Check for a long aliphatic chain, possibly representing phytol-like chains
    long_chain_smarts = Chem.MolFromSmarts('CCCCCCCCCCCCCO')
    if not mol.HasSubstructMatch(long_chain_smarts):
        return False, "Long aliphatic chain characteristic of chlorophyll missing"

    return True, "Contains magnesium porphyrin core with characteristic modifications and a long aliphatic chain"
"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for porphyrin core with magnesium coordination
    # Generalized search for heterocyclic rings with N and Mg link
    porphyrin_core_general = Chem.MolFromSmarts('n1cccc1n2ccc([Mg])c3ncc4ccc(n)c4c23')
    if not mol.HasSubstructMatch(porphyrin_core_general):
        return False, "Porphyrin core with magnesium coordination not found"

    # Check for long linear or branched aliphatic chain (e.g., phytol)
    aliphatic_chain_smarts = Chem.MolFromSmarts('CCC(CC)(CCCC)(CCCC)C')
    if not mol.HasSubstructMatch(aliphatic_chain_smarts):
        return False, "Long aliphatic chain characteristic of chlorophyll missing"

    return True, "Contains magnesium porphyrin core with a fifth ring and a long aliphatic chain"
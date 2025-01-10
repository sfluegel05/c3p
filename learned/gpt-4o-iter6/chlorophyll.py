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

    # Check for porphyrin core containing four pyrrole-like nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N')
    if n_count < 4:
        return False, "Less than four nitrogen atoms found"

    # Using SMARTS to look for typical chlorophyll porphyrin rings fused with the fifth ring
    porphyrin_core = Chem.MolFromSmarts('N1C=C2[Mg]53NC3=CC4=N[CH]C5=CC6=NC(=C2)C=1[C@@H](C[CHO]C(=O)O)C4=N6')
    if not mol.HasSubstructMatch(porphyrin_core):
        return False, "Typical chlorophyll porphyrin core with a fifth ring not found"

    # Check for long hydrophobic chains (e.g., phytol). This typically means many single-bond chain segments
    aliphatic_chain_smarts = Chem.MolFromSmarts('CCCCCCCCCCC')
    if not mol.HasSubstructMatch(aliphatic_chain_smarts):
        return False, "No long aliphatic chain characteristic of chlorophyll found"

    return True, "Contains magnesium porphyrin core with a fifth ring and a long aliphatic chain"
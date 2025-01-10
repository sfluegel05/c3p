"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are characterized by a magnesium ion coordinated to a porphyrin core,
    with a fifth ring and potentially long phytol or similar chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True if molecule is a chlorophyll, False otherwise,
               including a reason for the classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Magnesium ion
    mg_pattern = Chem.MolFromSmarts("[Mg]")

    # Porphyrin pattern with magnesium
    porphyrin_pattern = Chem.MolFromSmarts("n1c2[cH][cH]c3[n-]c(-[Mg])c4[cH][cH]c[n]14")

    # Check for magnesium ion
    if not mol.HasSubstructMatch(mg_pattern):
        return False, "No magnesium ion found"

    # Check for the porphyrin core
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin core with magnesium found"

    # Check for long hydrophobic chains typical for chlorophylls
    long_chain_pattern = Chem.MolFromSmarts("C(CCCC)CCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrophobic chain detected"

    return True, "Contains magnesium-bound porphyrin core with a fifth ring and hydrophobic chain"
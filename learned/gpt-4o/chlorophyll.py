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

    # Basic porphyrin pattern with magnesium - checking for tetrapyrrole with central Mg
    porphyrin_core_pattern = Chem.MolFromSmarts("[n]1[n][NH][n]c1-[Mg]-[n]2[n][NH][n]c2") 

    # Check for magnesium ion presence
    if not mol.HasSubstructMatch(mg_pattern):
        return False, "No magnesium ion found"

    # Check for the porphyrin core with magnesium
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
        return False, "No porphyrin core with magnesium found"

    # Check for the fifth ring beyond the porphyrin core
    fifth_ring_pattern = Chem.MolFromSmarts("c12cccc(nc1nccc2)") 
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring found in structure"

    # Check for long hydrophobic chains typical of chlorophylls
    phytol_chain_pattern = Chem.MolFromSmarts("C(=O)OCC(C)CC(C)CCCC(C)C")
    if not mol.HasSubstructMatch(phytol_chain_pattern):
        return False, "No long phytol chain or similar hydrophobic chain detected"
    
    return True, "Contains magnesium-bound porphyrin core with a fifth ring and hydrophobic chain"
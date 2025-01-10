"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are characterized by a magnesium ion coordinated to a porphyrin-type core
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

    # Check for the presence of magnesium
    mg_pattern = Chem.MolFromSmarts("[Mg]")
    if not mol.HasSubstructMatch(mg_pattern):
        return False, "No magnesium ion found"
    
    # Check for a generic porphyrin-like structure with magnesium (more generalized)
    porphyrin_core_pattern = Chem.MolFromSmarts("[n][c]1[c][c]2[c][c]3[n][c]4[c][c]5[n][c]6[n][c]7[c][c]8[n]9-[Mg]-[n]10c8c1[c]8c6c(c53)c(c427)9c8102")
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
        return False, "No porphyrin-like core with magnesium found"
    
    # Check for a fifth ring fused to the macrocycle
    # Looking for varying sizes of additional aromatic/heterocyclic rings
    fifth_ring_pattern = Chem.MolFromSmarts("c12[n][c][c](c[n][c1)cc2")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring found connected to the core"

    # Check for long hydrophobic chains (more generic than phytol)
    long_chain_pattern = Chem.MolFromSmarts("C[C@H](CCCCCCCC)C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrophobic chain detected"

    return True, "Contains magnesium-bound porphyrin-like core with a fifth ring and hydrophobic chain"
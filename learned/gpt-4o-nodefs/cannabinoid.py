"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a potential cannabinoid based on its SMILES string.
    Cannabinoids are known for having aromatic structures, long alkyl chains,
    and potentially complex ring systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define cannabinoid-associated structural motifs
    phenolic_ring = Chem.MolFromSmarts("c1cc(O)ccc1")
    long_alkyl_chain = Chem.MolFromSmarts("C(C)CCCCCCCC")
    cyclic_ethers = Chem.MolFromSmarts("c1occc1")
    basic_aromatic = Chem.MolFromSmarts("c1ccccc1")
    
    # Check for structural motifs commonly found in cannabinoids
    phenolic_match = mol.HasSubstructMatch(phenolic_ring)
    alkyl_chain_match = mol.HasSubstructMatch(long_alkyl_chain)
    ether_match = mol.HasSubstructMatch(cyclic_ethers)
    aromatic_match = mol.HasSubstructMatch(basic_aromatic)

    if (phenolic_match and alkyl_chain_match) or (ether_match and alkyl_chain_match) or aromatic_match:
        reason = "Potential cannabinoid structure detected with:"
        if phenolic_match:
            reason += " phenolic ring"
        if alkyl_chain_match:
            reason += " long alkyl chain"
        if ether_match:
            reason += " cyclic ether group"
        if aromatic_match:
            reason += " aromatic structure"
        return True, reason

    return False, "No specific cannabinoid structure detected"
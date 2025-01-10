"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by aromatic rings, long hydrocarbon chains, and oxygen-containing functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for aromatic systems (e.g., phenolic or indole rings)
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic system characteristic found"
    
    # Check for presence of oxygen atoms (indicating potential functional groups)
    oxygen_pattern = Chem.MolFromSmarts("[#8]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen-containing functional groups found"
    
    # Check for long aliphatic chains (substantially long CH-chains)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chains detected"
    
    return True, "Aromatic rings, long hydrocarbon chains, and oxygen functional groups are present"

# Example SMILES strings for testing
smiles_examples = [
    "C1(=C(C=C(CCCCC)C=C1O)O)C/C=C(/CCC=C(C)C)\C",  # cannabigerol
    "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole
]

for smiles in smiles_examples:
    is_cannabinoid_result, reason = is_cannabinoid(smiles)
    print(f"SMILES: {smiles}, Is cannabinoid? {is_cannabinoid_result}, Reason: {reason}")
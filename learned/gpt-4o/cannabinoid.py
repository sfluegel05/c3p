"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by aromatic rings with attached complex hydrocarbon chains, oxygen in heterocycles or functional groups.

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
    
    # Pattern for aromatic system (such as phenolic rings)
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic system characteristic found"
    
    # Check for oxygen-containing functional groups
    oxygen_pattern = Chem.MolFromSmarts("[#8]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen-containing functional groups found"
    
    # Check for long aliphatic chains typically found in cannabinoids
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain detected"

    # Flexible pattern to identify potential cannabinoids, accounting for structure diversity
    cannabinoid_flex_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]~[*]~[CX4,CX3]~[CX4,CX3]~[OX2H1,OX2,oxide]")  
    if mol.HasSubstructMatch(cannabinoid_flex_pattern):
        return True, "Contains aromatic systems with associated cannabinoid structures"

    return False, "Does not meet cannabinoid criteria"

# Example SMILES strings for testing
smiles_examples = [
    "C1(=C(C=C(CCCCC)C=C1O)O)C/C=C(/CCC=C(C)C)\C",  # cannabigerol
    "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole
]

for smiles in smiles_examples:
    is_cannabinoid_result, reason = is_cannabinoid(smiles)
    print(f"SMILES: {smiles}, Is cannabinoid? {is_cannabinoid_result}, Reason: {reason}")
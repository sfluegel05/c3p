"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by pharmacological activity and presence of oxygen in heterocyclic rings or as functional groups.

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

    # Check for long hydrocarbon chain patterns
    long_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C~C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for oxygen in functional groups or heterocyclic rings
    oxygen_pattern = Chem.MolFromSmarts("[#8]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygen-containing functional groups found"

    # Optional: Check for heterocycles (aromatic rings that may include oxygen)
    heterocycle_pattern = Chem.MolFromSmarts("[r6,r7,r8]")
    if not mol.HasSubstructMatch(heterocycle_pattern):
        return False, "No heterocyclic rings found"

    # Additional feature check for diverse cannabinoid structures
    # Here, consider adding more specific substructure searches based on known cannabinoids
    
    return True, "Contains characteristic long hydrocarbon chain and oxygen-containing functional groups"

# Example SMILES strings for testing
smiles_examples = [
    "C1(=C(C=C(CCCCC)C=C1O)O)C/C=C(/CCC=C(C)C)\C",  # cannabigerol
    "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole
]

for smiles in smiles_examples:
    is_cannabinoid_result, reason = is_cannabinoid(smiles)
    print(f"SMILES: {smiles}, Is cannabinoid? {is_cannabinoid_result}, Reason: {reason}")
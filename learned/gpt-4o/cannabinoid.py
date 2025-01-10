"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by complex hydrocarbon chains and oxygen within specific ring systems or as functional groups,
    typically linked to Cannabis plant metabolites and associated products.

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

    # Check improved long chain patterns with flexibility in structure
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Needs other contextual checks
    if mol.HasSubstructMatch(long_chain_pattern):
        # Check for oxygen functional groups
        oxygen_pattern = Chem.MolFromSmarts("[#8]")
        if not mol.HasSubstructMatch(oxygen_pattern):
            return False, "Detected long chain but no oxygen-containing functional groups found"

        # Check for at least one heterocycle 
        heterocycle_pattern = Chem.MolFromSmarts("[n,r;!R0]")
        if not mol.HasSubstructMatch(heterocycle_pattern):
            return False, "No heterocyclic rings or accepted substituents found"

        # Specific cannabinoid substructure patterns (potentially based on known variations)
        cannabinoid_substructure_pattern = Chem.MolFromSmarts("C1(CCCCC1)C=CC")
        if mol.HasSubstructMatch(cannabinoid_substructure_pattern):
            return True, "Contains characteristic cannabinoid substructures"

    return False, "Does not meet cannabinoid criteria"

# Example SMILES strings for testing
smiles_examples = [
    "C1(=C(C=C(CCCCC)C=C1O)O)C/C=C(/CCC=C(C)C)\C",  # cannabigerol
    "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole
]

for smiles in smiles_examples:
    is_cannabinoid_result, reason = is_cannabinoid(smiles)
    print(f"SMILES: {smiles}, Is cannabinoid? {is_cannabinoid_result}, Reason: {reason}")
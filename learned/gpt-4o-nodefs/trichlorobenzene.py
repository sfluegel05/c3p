"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene is defined as a benzene ring with three chlorine atoms substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for benzene ring with three chlorine atoms
    trichlorobenzene_pattern = Chem.MolFromSmarts("c1cc(Cl)cc(Cl)c(Cl)c1")
    
    if mol.HasSubstructMatch(trichlorobenzene_pattern):
        return True, "Contains a benzene ring with three chlorine atoms"

    return False, "Does not contain a trichlorobenzene structure"


# Example usage
smiles_strings = [
    "Clc1cc(Cl)cc(Cl)c1",     # 1,3,5-trichlorobenzene
    "Clc1cccc(Cl)c1Cl",       # 1,2,3-trichlorobenzene
    "Clc1ccc(Cl)c(Cl)c1",     # 1,2,4-trichlorobenzene
    "COc1cc(Cl)c(Cl)cc1Cl"    # Not strictly trichlorobenzene due to additional functional groups
]

for smiles in smiles_strings:
    result, reason = is_trichlorobenzene(smiles)
    print(f"SMILES: {smiles} => Result: {result}, Reason: {reason}")
"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene motif based on its SMILES string.
    A trichlorobenzene motif is defined as a benzene ring with exactly three chlorine atom substituents,
    potentially part of larger molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a trichlorobenzene motif, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define pattern for a benzene ring with exactly three chlorine atom substituents
    trichlorobenzene_pattern = Chem.MolFromSmarts("c1c([Cl])c([Cl])c([Cl])c(*)c1")
    if mol.HasSubstructMatch(trichlorobenzene_pattern):
        return True, "Contains a trichlorobenzene motif (benzene ring with three chlorine atoms)"

    return False, "Does not contain a trichlorobenzene motif"
    
# Example usage
smiles_strings = [
    "Clc1cc(Cl)cc(Cl)c1",     # 1,3,5-trichlorobenzene
    "Clc1cccc(Cl)c1Cl",       # 1,2,3-trichlorobenzene
    "Clc1ccc(Cl)c(Cl)c1",     # 1,2,4-trichlorobenzene
    "COc1cc(Cl)c(Cl)cc1Cl"    # 2,4,6-trichloroanisole
]

for smiles in smiles_strings:
    result, reason = is_trichlorobenzene(smiles)
    print(f"SMILES: {smiles} => Result: {result}, Reason: {reason}")
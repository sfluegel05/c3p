"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:52217 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hemiaminal pattern: a carbon atom bonded to exactly one amino group and one hydroxy group
    # The carbon can be part of a ring or have double bonds
    hemiaminal_pattern = Chem.MolFromSmarts("[C;H0]([N;H1,H2])([O;H1])")
    
    # Check if the molecule contains the hemiaminal pattern
    if mol.HasSubstructMatch(hemiaminal_pattern):
        # Further check to ensure the carbon is not part of a more stable functional group
        # Exclude amides, esters, carbamates, and other stable groups
        amide_pattern = Chem.MolFromSmarts("[C;H0]([N;H1,H2])([O;H1])([C;H0]=O)")
        ester_pattern = Chem.MolFromSmarts("[C;H0]([O;H1])([C;H0]=O)")
        carbamate_pattern = Chem.MolFromSmarts("[C;H0]([N;H1,H2])([O;H1])([C;H0]=O)")
        
        if not mol.HasSubstructMatch(amide_pattern) and not mol.HasSubstructMatch(ester_pattern) and not mol.HasSubstructMatch(carbamate_pattern):
            return True, "Contains a carbon atom bonded to both an amino group and a hydroxy group"
        else:
            return False, "Carbon atom is part of a more stable functional group (e.g., amide, ester, carbamate)"
    else:
        return False, "No carbon atom found with both an amino group and a hydroxy group attached"
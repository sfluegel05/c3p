"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is characterized by the -N=C=S group, where the N is directly attached to a carbon or heteroatom (not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for isothiocyanate group, ensuring N is bonded to a non-H atom or a heteroatom
    isothiocyanate_pattern = Chem.MolFromSmarts("[!H;!C;!N;!O;!S;!P][N]=[C]=[S]")
    
    # Check for substructure matches
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)

    if len(matches) == 1 : #Expect only one match for this class
        return True, "Contains isothiocyanate group"
    elif len(matches) > 1 :
        return False, f"Contains more than one isothiocyanate group: {len(matches)} matches"
    else:
        return False, "Missing isothiocyanate group"
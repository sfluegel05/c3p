"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:23003 chalcone
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core chalcone pattern: Ar-CH=CH-C(=O)-Ar
    # This pattern allows for substitutions and heteroatoms on the aromatic rings
    chalcone_pattern = Chem.MolFromSmarts("[*]1[*][*][*][*]1-[CH]=[CH]-[C](=O)-[*]2[*][*][*][*]2")
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "Core chalcone structure (Ar-CH=CH-C(=O)-Ar) not found"

    # Ensure that the aromatic rings are properly connected to the propenone structure
    matches = mol.GetSubstructMatches(chalcone_pattern)
    for match in matches:
        # Get the atoms in the match
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Check that the aromatic rings are connected to the propenone structure
        if (atoms[0].GetSymbol() in ['C', 'N', 'O'] and atoms[5].GetSymbol() in ['C', 'N', 'O'] and
            atoms[0].GetDegree() >= 2 and atoms[5].GetDegree() >= 2):
            return True, "Contains core chalcone structure (Ar-CH=CH-C(=O)-Ar) with possible substitutions"

    return False, "Aromatic rings not properly connected to the propenone structure"
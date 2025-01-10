"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of the -NH2 group attached to a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for primary amine pattern: a nitrogen with exactly two hydrogens
    # attached to a single carbon (hydrocarbyl group)
    primary_amine_pattern = Chem.MolFromSmarts("[NH2][CX4]")  
    # NH2: nitrogen with exactly two hydrogens
    # CX4: carbon with four single bonds (aliphatic)

    matches = mol.GetSubstructMatches(primary_amine_pattern)
    
    # Ensure the NH2 group is not linked to more complex structures indicative of non-primary amines
    if matches:
        # Verify not bonded to additional non-carbon atoms that may define a secondary/tertiary amine 
        for match in matches:
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            if not any(b.GetOtherAtomIdx(nitrogen_atom.GetIdx()) != match[1] 
                       for b in nitrogen_atom.GetBonds()):
                return True, "Contains primary amine group (-NH2) attached to a hydrocarbyl group"

    return False, "Primary amine group not found"
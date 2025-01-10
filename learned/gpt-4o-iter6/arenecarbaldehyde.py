"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as any aldehyde with the carbonyl group attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined aldehyde pattern where the carbonyl carbon is attached to a benzene-like aromatic carbon
    aldehyde_aromatic_attachment = Chem.MolFromSmarts("[C](=[O])[c;R1]")

    if mol.HasSubstructMatch(aldehyde_aromatic_attachment):
        # Additional check to ensure the aromatic nature of the aldehyde attachment
        aromatic_atom_ids = [m[1] for m in mol.GetSubstructMatches(aldehyde_aromatic_attachment)]
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in aromatic_atom_ids):
            return True, "Aldehyde group attached to an aromatic moiety"
    
    return False, "Aldehyde group not attached to aromatic moiety"
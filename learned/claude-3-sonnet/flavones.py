"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:28167 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is a flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-arylchromen-4-one skeleton
    # 1. Find the chromone core
    chromone_pattern = Chem.MolFromSmarts("O=C1C=C(Oc2ccccc2)C=CC1")
    chromone_matches = mol.GetSubstructMatches(chromone_pattern)
    if not chromone_matches:
        return False, "No chromone core found"
    
    # 2. Check for aromatic ring attached to the chromone
    for match in chromone_matches:
        chromone_atom = mol.GetAtomWithIdx(match[2])
        for neighbor in chromone_atom.GetNeighbors():
            ring_info = mol.GetAtomRingInfo(neighbor)
            for ring in ring_info.AtomRings():
                if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True, "Contains 2-arylchromen-4-one skeleton"
    
    return False, "No 2-arylchromen-4-one skeleton found"
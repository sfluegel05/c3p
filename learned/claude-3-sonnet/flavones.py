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

    # Look for chromone core
    chromone_pattern = Chem.MolFromSmarts("C1=CC(=O)Oc2ccccc12")
    chromone_matches = mol.GetSubstructMatches(chromone_pattern)
    if not chromone_matches:
        return False, "No chromone core found"
    
    # Look for phenyl ring attached to the chromone
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
    if not phenyl_matches:
        return False, "No phenyl ring found"
    
    # Check if the phenyl ring is attached to the chromone core
    ring_info = mol.GetRingInfo()
    for chromone_match in chromone_matches:
        chromone_atom = mol.GetAtomWithIdx(chromone_match[2])
        for neighbor in chromone_atom.GetNeighbors():
            if neighbor.GetIsAromatic() and any(ring_info.IsAtomInRingOfSize(neighbor_idx, 6) for neighbor_idx in neighbor.GetNeighbors()):
                break
        else:
            continue
        break
    else:
        return False, "No 2-arylchromen-4-one skeleton found"
    
    # Check for common substituents
    has_hydroxyl = any(atom.GetSymbol() == "O" and atom.GetHybridization() == Chem.HybridizationType.SP3 for atom in mol.GetAtoms())
    has_methoxy = any(atom.GetSymbol() == "O" and any(neighbor.GetSymbol() == "C" for neighbor in atom.GetNeighbors()) for atom in mol.GetAtoms())
    
    if has_hydroxyl or has_methoxy:
        return True, "Contains 2-arylchromen-4-one skeleton with hydroxyl or methoxy substituents"
    else:
        return True, "Contains 2-arylchromen-4-one skeleton"
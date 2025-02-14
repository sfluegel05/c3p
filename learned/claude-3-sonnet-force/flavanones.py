"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:26176 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone is a flavan with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for core flavanone scaffold
    scaffold_pattern = Chem.MolFromSmarts("[OX2]C1=C(C=C2C=CC=CC=2)C(=O)CCC1")
    scaffold_match = mol.GetSubstructMatch(scaffold_pattern)
    if not scaffold_match:
        return False, "No valid flavanone scaffold found"

    # Check for pyran ring
    pyran_ring = mol.GetSubstructMatch(Chem.MolFromSmarts("[OX2]C1=CC(=O)CCC1"))
    if not pyran_ring:
        return False, "No pyran ring found"

    # Check for benzene ring
    benzene_ring = mol.GetSubstructMatch(Chem.MolFromSmarts("c1ccccc1"))
    if not benzene_ring:
        return False, "No benzene ring found"

    # Check for connectivity between pyran and benzene rings
    conn_pattern = Chem.MolFromSmarts("[OX2]C1=CC(=O)CCC1C2=CC=CC=C2")
    conn_match = mol.GetSubstructMatch(conn_pattern)
    if not conn_match:
        return False, "Pyran and benzene rings not connected"

    # Check for common substituents or structural features
    has_hydroxy_groups = any(atom.GetAtomicNum() == 8 and atom.GetIsAromatic() is False for atom in mol.GetAtoms())
    has_methoxy_groups = any(atom.GetAtomicNum() == 8 and atom.GetIsAromatic() is False and any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()) for atom in mol.GetAtoms())
    has_prenyl_groups = any(atom.GetAtomicNum() == 6 and atom.GetDegree() == 4 and any(neighbor.GetAtomicNum() == 8 for neighbor in atom.GetNeighbors()) and any(neighbor.GetAtomicNum() == 6 and len(Chem.FindAllPathsOfLengthN(mol, neighbor.GetIdx(), 3, useBondOrder=False)) for neighbor in atom.GetNeighbors()) for atom in mol.GetAtoms())

    if not has_hydroxy_groups and not has_methoxy_groups and not has_prenyl_groups:
        return False, "No common substituents or structural features found"

    return True, "Contains flavanone scaffold with pyran and benzene rings, and common substituents or structural features"
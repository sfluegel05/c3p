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

    # Check for benzene ring
    benzene_ring = mol.GetSubstructMatch(Chem.MolFromSmarts("c1ccccc1"))
    if not benzene_ring:
        return False, "No benzene ring found"

    # Check for pyran ring with carbonyl and ether groups
    pyran_pattern = Chem.MolFromSmarts("[OX2]C1=C(O)C(=O)CCC1")
    pyran_match = mol.GetSubstructMatch(pyran_pattern)
    if not pyran_match:
        return False, "No pyran ring with carbonyl and ether groups found"

    # Check for connectivity between benzene and pyran rings
    flavanone_scaffold = mol.GetSubstructMatch(Chem.MolFromSmarts("[OX2]C1=C(O)C(=O)CCC1C2=CC=CC=C2"))
    if not flavanone_scaffold:
        return False, "No valid flavanone scaffold found"

    # Additional checks for common substituents or structural features (optional)
    # ...

    return True, "Contains flavanone scaffold with benzene and pyran rings, carbonyl and ether groups"
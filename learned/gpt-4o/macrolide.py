"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Macrocyclic lactone pattern: detect a large ring containing the ester group
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)CCCCCCCCCCC1")
    
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No large macrocyclic lactone group found"

    # Check cycle sizes, ensure the presence of a macrocycle with the pattern
    ssr = Chem.GetSymmSSSR(mol)
    large_rings = [ring for ring in ssr if len(ring) >= 12]
    
    for ring in large_rings:
        submol = Chem.PathToSubmol(mol, ring)
        if submol.HasSubstructMatch(lactone_pattern):
            return True, f"Contains macrocyclic lactone ring of size {len(ring)}"
    
    return False, "No sufficiently large macrocyclic lactone group found"
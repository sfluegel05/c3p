"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent attached to the phenol ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible phenol group pattern to account for potential substituents
    phenol_pattern = Chem.MolFromSmarts("c1ccc([#6,O,N;R0,R1,R2,R3;!$([OH])]c1)O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) < 1:
        return False, "No phenol group found"

    # Find nitro groups attached to aromatic carbon
    nitro_pattern = Chem.MolFromSmarts("[$(c[N+](=O)[O-])]")
    nitro_matches = [match for match in mol.GetSubstructMatches(nitro_pattern)]
    
    # Check if any nitro group is directly attached to the detected phenol group(s)
    if not any(any(mol.GetAtomWithIdx(n[0]).GetIdx() in set(sum(phenol_matches, ())) for n in match) for match in nitro_matches):
        return False, f"No nitro group directly attached to phenol group"

    return True, "Contains phenol group with a nitro group directly attached"
"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is chlorophyll based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for magnesium (Mg) ion coordinated in the molecule
    mg_atom = any(atom.GetAtomicNum() == 12 for atom in mol.GetAtoms())
    if not mg_atom:
        return False, "No magnesium atom found"

    # Check for a simplified macrocycle structure (porphyrin/chlorin-like)
    # Chlorin is a more reduced form of porphyrin, so a simpler pattern
    chlorin_pattern = Chem.MolFromSmarts("c1cc[nH]c2ccc(nc3cc[nH]c(cc4[nH]c5ccccc5[nH]4)c23)c1")
    if not mol.HasSubstructMatch(chlorin_pattern):
        return False, "No chlorin-like macrocyclic structure found"

    # Check for presence of hydrocarbon chains that resemble the phytol group
    # More generalized pattern to match varied hydrocarbon tails
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C=C([C@H](CC)C)C")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No suitable hydrocarbon chains found"

    return True, "Matches chlorophyll characteristics: macrocyclic structure with magnesium and hydrocarbon chains"
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

    # Check for chlorin-like structure (porphyrin with reduced bonds in one ring)
    chlorin_pattern = Chem.MolFromSmarts("n1c(c2nc(nc2c3cc1[nH+]3)C)c4[nH]c5c(n4)[nH+]5")
    if not mol.HasSubstructMatch(chlorin_pattern):
        return False, "No chlorin-like structure found"

    # Check for long hydrocarbon chains or common hydrocarbon fragments in chlorophyll
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No suitable long hydrocarbon chains found"
    
    # Check for presence of ester groups
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found, common in chlorophyll"

    return True, "Matches chlorophyll characteristics: chlorin-like structure with magnesium, hydrocarbon chains, and ester groups"
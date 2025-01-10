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

    # Improved pattern for detecting chlorin-like macrocyclic structure
    chlorin_pattern = Chem.MolFromSmarts("n1c(c2nc(cc3c2c1)c(c4cc(c(c5[nH]c(nc5c6ccc(n6)c3)c4))c7[nH]c8c(c7)c(c9cc(nc(n9)c8)c0ccc(nc0)c1)c2)n3)cc4)c2n8c(c)c(=n9)ccn(c9)cc8c3c(c4)n(cc3)c2")
    if not mol.HasSubstructMatch(chlorin_pattern):
        return False, "No chlorin-like macrocyclic structure found"

    # Check for presence of phytol or geranylgeranyl-like hydrocarbon chains commonly found
    phytol_chain_pattern = Chem.MolFromSmarts("C(CC)=C/C\\C=C/C\\C=C/C")
    if not mol.HasSubstructMatch(phytol_chain_pattern):
        return False, "No suitable phytol or similar hydrocarbon chains found"

    # Expanded pattern for presence of ester groups in given positions
    ester_pattern = Chem.MolFromSmarts("CO[C](=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups in appropriate locations found"

    return True, "Matches chlorophyll characteristics: macrocyclic structure with magnesium, phytol-like chains, and ester groups"
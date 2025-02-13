"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for at least one hydroxy group (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy (OH) group found"
    
    # Check the carbon chain length; lower threshold to include smaller fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:  # Reduce minimum carbon count for hydroxy fatty acids
        return False, "Too few carbon atoms for a fatty acid"
    
    # Allow limited branching or rings, but ensure molecule is primarily a linear fatty acid
    ring_info = mol.GetRingInfo()
    max_rings = 2  # Allow some structures like cyclic ethers
    if ring_info.NumRings() > max_rings:
        return False, f"More than {max_rings} rings detected, inconsistent with typical fatty acid structure"
    
    # Check for peptide-like structures as a further refinement
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Peptide-like structures detected"

    return True, "Molecule contains a carboxylic acid group and one or more hydroxy groups; matches hydroxy fatty acid structure"
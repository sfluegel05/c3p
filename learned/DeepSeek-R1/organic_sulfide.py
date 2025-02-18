"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:48366 organic sulfide (RSR', R ≠ H)
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    Organic sulfides have the structure R-S-R' where R and R' are non-hydrogen atoms connected via single bonds,
    and sulfur is not part of an aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            # Skip aromatic sulfurs (e.g., in thiophene)
            if atom.GetIsAromatic():
                continue
                
            # Check sulfur has exactly two single bonds
            if atom.GetDegree() != 2:
                continue
                
            bonds = atom.GetBonds()
            # Verify all bonds are single bonds (no double bonds or aromatic)
            if any(bond.GetBondType() != Chem.BondType.SINGLE for bond in bonds):
                continue
                
            # Check both neighbors are non-hydrogen
            neighbors = atom.GetNeighbors()
            if all(n.GetAtomicNum() != 1 for n in neighbors):
                return True, "Contains sulfide group (R-S-R') with both R groups non-hydrogen via single bonds"
    
    return False, "No sulfide group (R-S-R') with both R groups non-hydrogen found"
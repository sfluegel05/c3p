"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule contains a tertiary amine group.
    A tertiary amine has a nitrogen atom with three carbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Skip if nitrogen is positively charged
            if atom.GetFormalCharge() > 0:
                continue
                
            # Skip if nitrogen is in a ring and has aromatic bonds
            if atom.IsInRing() and atom.GetIsAromatic():
                continue
                
            # Count carbon neighbors
            carbon_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    carbon_neighbors += 1
            
            # Check if it's connected to exactly 3 carbons
            if carbon_neighbors == 3:
                # Check if it's not part of an amide group
                is_amide = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        # Check if this carbon is part of C=O
                        for sub_neighbor in neighbor.GetNeighbors():
                            if sub_neighbor.GetAtomicNum() == 8 and sub_neighbor.GetBonds()[0].GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                is_amide = True
                                break
                
                if not is_amide:
                    return True, "Contains nitrogen atom bonded to exactly three carbon atoms"

    return False, "No tertiary amine group found"
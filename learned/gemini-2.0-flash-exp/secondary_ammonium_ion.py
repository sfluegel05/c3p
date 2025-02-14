"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is a nitrogen atom with a positive charge and two directly attached
    carbon containing groups via single bonds, and not part of any ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1: # Check for nitrogen and positive charge
            if atom.IsInRing(): # Exclude nitrogens in rings
                continue
            carbon_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == rdchem.BondType.SINGLE: # Check if neighbour is a carbon via single bond
                    carbon_neighbors += 1
            if carbon_neighbors == 2:
                return True, "Contains at least one positively charged nitrogen with two directly attached carbon groups via single bonds, and not part of any ring."


    return False, "Does not contain a positively charged nitrogen with two directly attached carbon groups via single bonds, or present in ring."
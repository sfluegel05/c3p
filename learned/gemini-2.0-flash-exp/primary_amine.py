"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is defined as a nitrogen atom bonded to one other non-hydrogen substituent and 0, 1, or 2 hydrogen atoms (implied or explicit).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES, remove any salt info using RDKit fragment extraction
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get largest fragment if multiple fragments are present
    frags = rdmolops.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Check each nitrogen atom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7: # Check if the atom is Nitrogen
            num_non_h = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() != 1: # If neighbor is not H
                    num_non_h += 1

            if num_non_h == 1 and atom.GetFormalCharge() == 0 and not atom.IsInRing():
                  return True, "Contains a primary amine"

    return False, "Does not contain a primary amine"
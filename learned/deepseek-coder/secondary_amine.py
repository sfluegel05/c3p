"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32877 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary amine
    # The pattern looks for a nitrogen atom bonded to exactly two carbon atoms and one hydrogen atom
    # and ensures that the nitrogen is not part of a double bond or a heterocycle unless it is a secondary amine
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])([#6])")

    # Check if the molecule contains the secondary amine pattern
    if mol.HasSubstructMatch(secondary_amine_pattern):
        # Further check to exclude amides, imines, and other non-secondary amine contexts
        for match in mol.GetSubstructMatches(secondary_amine_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            # Ensure the nitrogen is not part of a double bond
            if not any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in nitrogen_atom.GetBonds()):
                # Ensure the nitrogen is not part of a heterocycle unless it is a secondary amine
                if not nitrogen_atom.IsInRing() or (nitrogen_atom.IsInRing() and nitrogen_atom.GetDegree() == 2):
                    return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom (secondary amine)"
    
    return False, "No secondary amine pattern found"
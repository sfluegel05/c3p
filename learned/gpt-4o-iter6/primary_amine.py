"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of an -NH2 group attached to a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific SMARTS pattern for discriminating primary amines
    # Primary amines are NH2 groups attached to carbons, ensuring hydrocarbyl chains allowed, avoid adjacent heteroatoms
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][CX3,CX4,c]")

    matches = mol.GetSubstructMatches(primary_amine_pattern)

    # Analyze matches for presence of specific attachment criteria
    for match in matches:
        nitrogen_atom = mol.GetAtomWithIdx(match[0])
        bonded_carbons = [b.GetOtherAtom(nitrogen_atom) for b in nitrogen_atom.GetBonds() if b.GetOtherAtom(nitrogen_atom).GetAtomicNum() == 6]

        for carbon in bonded_carbons:
            # Ensure that nitrogen is bonded to a carbon which has valid hydrocarbyl connections
            # This ensures primary amine is part of a straightforward hydrocarbyl group
            hydrocarbon_bonds = any(b.GetOtherAtom(carbon).GetAtomicNum() == 6 for b in carbon.GetBonds() if b.GetOtherAtom(carbon) != nitrogen_atom)
            
            if hydrocarbon_bonds:
                return True, "Primary amine group (-NH2) correctly identified with hydrocarbyl connection"

    return False, "Primary amine group not found"
"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime contains an aliphatic chain and an aldehyde group
    converted to an oxime (R-C=NOH), without presence in aromatic systems or rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the oxime group pattern "-C=NO"
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OH]")
    oxime_matches = mol.GetSubstructMatches(oxime_pattern)
    if not oxime_matches:
        return False, "No oxime group found"

    # Ensure oxime is part of an aliphatic environment
    for match in oxime_matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        
        # Ensure the oxime carbon is sp3 hybridized and not in a ring
        if carbon_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or carbon_atom.IsInRing():
            return False, "Oxime group is part of an aromatic ring or highly conjugated system"

        # Check neighboring atoms to ensure they are part of an aliphatic chain
        is_aliphatic = True
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic() or neighbor.IsInRing():
                is_aliphatic = False
                break
        
        if is_aliphatic:
            return True, "Contains an aliphatic aldoxime group"

    return False, "Oxime group is not part of a purely aliphatic chain"
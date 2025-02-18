"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:39194 methyl sulfide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where at least one organyl group attached to sulfur is methyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for methyl sulfide where sulfur is connected to methyl and another aliphatic carbon
    # [SX2;!a] ensures non-aromatic thioether, [CH3;!a] ensures methyl isn't part of aromatic system
    sulfide_pattern = MolFromSmarts('[SX2;!a]([CH3;!a])[#6;!a]')
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Aliphatic sulfide with methyl group attached to sulfur"

    # Check for dimethyl sulfide case with both groups non-aromatic
    dimethyl_pattern = MolFromSmarts('[SX2;!a]([CH3;!a])([CH3;!a])')
    if mol.HasSubstructMatch(dimethyl_pattern):
        return True, "Dimethyl sulfide (both groups are methyl)"

    # Additional check for branched methyl groups or complex aliphatic chains
    # Iterate through all sulfur atoms to catch edge cases
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 16 and 
            not atom.GetIsAromatic() and 
            atom.GetDegree() == 2):
            
            neighbors = atom.GetNeighbors()
            methyl_found = any(
                n.GetAtomicNum() == 6 and 
                not n.GetIsAromatic() and 
                n.GetTotalNumHs() >= 3  # Indicates CH3 group
                for n in neighbors
            )
            
            if methyl_found:
                # Ensure other substituent is aliphatic carbon
                other_substituent = next(n for n in neighbors if not (n.GetAtomicNum() == 6 and not n.GetIsAromatic() and n.GetTotalNumHs() >= 3))
                if other_substituent.GetAtomicNum() == 6 and not other_substituent.GetIsAromatic():
                    return True, "Aliphatic sulfide with methyl group attached to sulfur"

    return False, "No aliphatic sulfide with methyl group found"
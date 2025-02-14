"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: CHEBI:24428 alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2,
    and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if SMILES contains only valid alkane atoms and bonds
    allowed_atoms = [6, 1]  # C and H
    allowed_bonds = [Chem.BondType.SINGLE]
    if any(atom.GetAtomicNum() not in allowed_atoms for atom in mol.GetAtoms()) or \
       any(bond.GetBondType() not in allowed_bonds for bond in mol.GetBonds()):
        return False, "Contains invalid atoms or bonds for an alkane"

    # Check for cycles
    if any(ring for ring in mol.GetRingInfo().AtomRings()):
        return False, "Contains cyclic structures"

    # Check for substituents and their connectivity
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            substituents = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if len(substituents) > 4:
                return False, "Carbon atom has more than 4 substituents"
            elif 6 in substituents and any(substituent != 1 for substituent in substituents):
                return False, "Contains non-alkyl substituents"

    return True, "Acyclic compound containing only carbon and hydrogen in the correct ratio"
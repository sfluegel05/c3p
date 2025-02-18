"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose (5 carbons) with a (potential) aldehyde group.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for exactly 5 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 5:
        return False, f"Expected 5 carbons, found {carbon_count}"

    # Check for explicit aldehyde group (O=CH-)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains aldehyde group"

    # Check for cyclic hemiacetal form (pentose sugar)
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        # Check if ring contains an oxygen
        has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
        if has_oxygen:
            # Count hydroxyl groups attached to ring carbons
            hydroxyl_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon in ring
                    # Check for adjacent hydroxyl groups
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                            hydroxyl_count += 1
                            break  # Count once per carbon
            if hydroxyl_count >= 3:
                return True, "Cyclic pentose with potential aldehyde group"

    return False, "Not an aldopentose"
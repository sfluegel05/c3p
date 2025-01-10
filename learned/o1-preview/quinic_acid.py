"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    Quinic acid is characterized by a cyclohexane ring with multiple hydroxyl or ester groups,
    and a carboxylic acid or ester group attached to the ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is quinic acid or its derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find cyclohexane rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    cyclohexane_rings = [ring for ring in rings if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]
    
    if not cyclohexane_rings:
        return False, "No cyclohexane ring found"

    # Check each cyclohexane ring
    for ring in cyclohexane_rings:
        oxy_substituents = 0
        carboxylic_acid_or_ester = False

        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetIdx() not in ring:
                    if nbr.GetAtomicNum() == 8:  # Oxygen
                        # Check if oxygen is double-bonded (carbonyl) or single-bonded (hydroxyl/ether)
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            oxy_substituents += 1
                        elif bond.GetBondType() == Chem.BondType.DOUBLE:
                            # Check for carboxylic acid or ester group
                            nbr2 = None
                            for b in nbr.GetBonds():
                                if b.GetOtherAtom(nbr).GetIdx() != atom.GetIdx():
                                    nbr2 = b.GetOtherAtom(nbr)
                                    break
                            if nbr2 and nbr2.GetAtomicNum() == 8:
                                carboxylic_acid_or_ester = True

        if oxy_substituents >= 3 and carboxylic_acid_or_ester:
            return True, "Contains quinic acid core structure"

    return False, "Does not match quinic acid core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'quinic acid',
        'definition': 'A cyclitol carboxylic acid.',
    },
}
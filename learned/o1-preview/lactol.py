"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.
    They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved SMARTS pattern for lactol (cyclic hemiacetal)
    # Pattern explanation:
    # - Ring structure with an ether oxygen: [O;R]
    # - Connected to ring carbon: [C;R]
    # - That carbon has a hydroxyl group: [C;R]([OH])
    # - This pattern allows for unsaturation and aromatic rings
    lactol_pattern = Chem.MolFromSmarts("[O;R][C;R]([OH])")
    
    if lactol_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for lactol substructure
    matches = mol.GetSubstructMatches(lactol_pattern)
    if matches:
        # Additional checks to exclude sugars and polysaccharides
        # Count the number of hydroxyl groups attached to ring carbons
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings > 0:
            for ring in ring_info.AtomRings():
                hydroxyls_on_ring = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    # Check if atom is a carbon with hydroxyl group
                    if atom.GetAtomicNum() == 6:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                                hydroxyls_on_ring +=1
                # If more than one hydroxyl on ring, likely a sugar, exclude
                if hydroxyls_on_ring >1:
                    return False, "Ring has multiple hydroxyl groups, likely a sugar"
        return True, "Contains lactol moiety (cyclic hemiacetal)"
    else:
        return False, "No lactol moiety found"

__metadata__ = {
    'chemical_class': {
        'name': 'lactol',
        'definition': 'Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group. They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.'
    },
}
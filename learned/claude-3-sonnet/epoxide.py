"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:36975 epoxide
An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an epoxide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for epoxide pattern
    epoxide_pattern = Chem.MolFromSmarts("[O;r3]")
    matches = mol.GetSubstructMatches(epoxide_pattern)

    if not matches:
        return False, "No epoxide group found"

    # Additional checks for valid epoxide
    for match in matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)

        # Check if oxygen atom is part of a 3-membered ring
        rings = mol.GetRingInfo().AtomRings()
        for ring in rings:
            if atom_idx in ring and len(ring) == 3:
                # Check bond orders and atom hybridization
                neighbors = [mol.GetAtomWithIdx(nb) for nb in atom.GetNeighbors()]
                if all(nb.GetHybridization() == Chem.HybridizationType.SP3 for nb in neighbors) and \
                   all(mol.GetBondBetweenAtoms(atom_idx, nb.GetIdx()).GetBondType() == Chem.BondType.SINGLE for nb in neighbors):
                    return True, "Molecule contains an epoxide group"

    return False, "No valid epoxide group found"
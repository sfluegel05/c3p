"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16646 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains exactly 6 carbon atoms and 6 oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count != 6:
        return False, "Not a hexose (does not contain 6 carbons and 6 oxygens)"

    # Check for the presence of a ring structure
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Check if there is a single ring of size 5 or 6 (pyranose or furanose)
    if len(rings) == 1 and (len(rings[0]) == 5 or len(rings[0]) == 6):
        # Identify the anomeric carbon (the carbon connected to the ring oxygen)
        anomeric_carbon = None
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 4:
                neighbors = [mol.GetAtomWithIdx(n).GetAtomicNum() for n in atom.GetNeighbors()]
                if 8 in neighbors:
                    anomeric_carbon = atom.GetIdx()
                    break

        if anomeric_carbon is not None:
            # Get the chiral centers of the molecule
            chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('_ChiralityPossible') and atom.GetProp('_ChiralityPossible') == 'Unk']

            # Check if there are at least 4 chiral centers (minimum for a hexose)
            if len(chiral_centers) >= 4:
                # Find the chiral center corresponding to position 5
                position_5 = None
                for idx in chiral_centers:
                    atom = mol.GetAtomWithIdx(idx)
                    ring_bond_count = sum(1 for neighbor in atom.GetNeighbors() if mol.GetAtomWithIdx(neighbor).IsInRing())
                    if ring_bond_count == 1:
                        position_5 = idx
                        break

                if position_5 is not None:
                    # Check if the chiral center at position 5 has D-configuration
                    atom = mol.GetAtomWithIdx(position_5)
                    if atom.GetProp('_CIPCode') == 'R':
                        return True, "Molecule has D-configuration at position 5"
                    else:
                        return False, "Molecule does not have D-configuration at position 5"

    # If no ring structure or linear form
    return False, "Not a hexose (no valid ring structure or linear form found)"
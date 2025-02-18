"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate has a six-membered carbon ring with myo-configuration (5 hydroxyls in one configuration,
    one axial) and at least one phosphate group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all six-membered carbon rings
    rings = mol.GetRingInfo().AtomRings()
    six_membered_carbon_rings = [ring for ring in rings if len(ring) == 6 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring)]

    if not six_membered_carbon_rings:
        return False, "No six-membered carbon ring"

    for ring in six_membered_carbon_rings:
        # Check all substituents on ring carbons are oxygen atoms
        all_oxygen_substituents = True
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring and neighbor.GetAtomicNum() != 8:
                    all_oxygen_substituents = False
                    break
            if not all_oxygen_substituents:
                break
        if not all_oxygen_substituents:
            continue

        # Check all six carbons are stereocenters
        stereocenters = []
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
                break
            stereocenters.append(atom.GetChiralTag())
        else:
            # Count CW and CCW
            count_cw = sum(1 for tag in stereocenters if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            count_ccw = sum(1 for tag in stereocenters if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

            if (count_cw == 5 and count_ccw == 1) or (count_ccw == 5 and count_cw == 1):
                # Check for at least one phosphate group (P connected via O)
                phosphate_found = False
                for a in ring:
                    atom = mol.GetAtomWithIdx(a)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            # Traverse from oxygen to find phosphorus
                            stack = [neighbor]
                            visited = set()
                            while stack and not phosphate_found:
                                current = stack.pop()
                                if current.GetIdx() in visited:
                                    continue
                                visited.add(current.GetIdx())
                                if current.GetAtomicNum() == 15:  # Phosphorus found
                                    phosphate_found = True
                                    break
                                # Add neighboring atoms except the ring carbon
                                for nbr in current.GetNeighbors():
                                    if nbr.GetIdx() != atom.GetIdx():
                                        stack.append(nbr)
                    if phosphate_found:
                        break
                if phosphate_found:
                    return True, "myo-inositol core with phosphate group(s)"

    return False, "Does not meet myo-inositol phosphate criteria"
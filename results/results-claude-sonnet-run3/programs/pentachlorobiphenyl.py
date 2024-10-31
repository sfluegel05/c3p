from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pentachlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a pentachlorobiphenyl (PCB with molecular formula C12H5Cl5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pentachlorobiphenyl, False otherwise 
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular formula by counting atoms
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != "C12H5Cl5":
        return False, f"Molecular formula must be C12H5Cl5, found {formula}"

    # Check for presence of two benzene rings
    rings = mol.GetRingInfo()
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if len(six_membered_rings) != 2:
        return False, "Must contain exactly 2 benzene rings"

    # Check for biphenyl structure (rings connected by single bond)
    ring1_atoms = set(six_membered_rings[0])
    ring2_atoms = set(six_membered_rings[1])
    
    # Find connecting atoms between rings
    connections = 0
    for atom1_idx in ring1_atoms:
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        for neighbor in atom1.GetNeighbors():
            if neighbor.GetIdx() in ring2_atoms:
                connections += 1

    if connections != 1:
        return False, "Rings must be connected by exactly one single bond"

    # All checks passed
    return True, "Valid pentachlorobiphenyl structure"
# Pr=1.0
# Recall=1.0
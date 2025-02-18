"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:46617 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of D-ribose connected to a nucleobase via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define flexible D-ribose pattern allowing substituents on hydroxyl positions
    ribose_pattern = Chem.MolFromSmarts(
        "[C@H]1(O[C@H](C[OX2])[C@H]([!C])[C@@H]([!C])1[!C])"
    )
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No D-ribose sugar detected"

    # Check for phosphate groups on the ribose's C5 (CH2OX group)
    for match in ribose_matches:
        c5_atom = match[2]  # Index of C in C[OX2] (CH2OX group)
        for neighbor in mol.GetAtomWithIdx(c5_atom).GetNeighbors():
            if neighbor.GetAtomicNum() == 15:  # Phosphorus atom
                return False, "Phosphate group present (nucleotide, not nucleoside)"

    # Check glycosidic bond (anomeric carbon connected to nucleobase heteroatom)
    glycosidic_bond = False
    nucleobase_atoms = set()
    for match in ribose_matches:
        anomeric_c = match[0]  # C1 in the ribose pattern
        for neighbor in mol.GetAtomWithIdx(anomeric_c).GetNeighbors():
            if neighbor.GetAtomicNum() in [7, 8]:  # N or O
                # Check if the heteroatom is part of a ring (crude nucleobase check)
                if neighbor.IsInRing():
                    glycosidic_bond = True
                    nucleobase_atoms.add(neighbor.GetIdx())
                    break
        if glycosidic_bond:
            break
    if not glycosidic_bond:
        return False, "No glycosidic bond to nucleobase heteroatom"

    # Check for common nucleobase patterns in the connected atoms
    nucleobase_patterns = [
        Chem.MolFromSmarts(sma) for sma in [
            "N1C=NC2=C1N=CN2",  # Adenine
            "N1C=NC2=C1NC(=O)N2",  # Guanine
            "C1=CNC(=O)NC1=O",  # Uracil
            "C1=C(NC(=O)N)C(=O)NC1",  # Thymine
            "N1C(=O)NC(=NC1)N"  # Cytosine
        ]
    ]
    has_base = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            has_base = True
            break
    if not has_base:
        # Check if any nucleobase atoms are aromatic (fallback)
        has_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in nucleobase_atoms)
        if not has_aromatic:
            return False, "No recognized nucleobase detected"

    return True, "D-ribose sugar connected to nucleobase via glycosidic bond"
"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate, obtained by
    esterification of glycerol with one fatty acid chain and phosphorylation on one hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define lysophosphatidic acid patterns
    # Acyl chain at sn-1 position, phosphate at sn-3
    lpa_sn1_pattern = Chem.MolFromSmarts("""
        [C:1]([O:4][P:5](=O)([O-])[O:6])([O:2][C:3])([H])
        (
            [O:2][C:3][O:7][C:8](=O)[C:9]
        )
    """)
    # Acyl chain at sn-2 position, phosphate at sn-3
    lpa_sn2_pattern = Chem.MolFromSmarts("""
        [C:1]([O:2][C:3])([O:4][P:5](=O)([O-])[O:6])([H])
        (
            [O:2][C:3][O:7][C:8](=O)[C:9]
        )
    """)
    # Acyl chain at sn-1 position, phosphate at sn-2
    lpa_sn1_phosphate_sn2 = Chem.MolFromSmarts("""
        [C:1]([O:2][P:5](=O)([O-])[O:6])([O:4][C:3])([H])
        (
            [O:4][C:3][O:7][C:8](=O)[C:9]
        )
    """)
    # Acyl chain at sn-2 position, phosphate at sn-1
    lpa_sn2_phosphate_sn1 = Chem.MolFromSmarts("""
        [C:1]([O:4][C:3])([O:2][P:5](=O)([O-])[O:6])([H])
        (
            [O:4][C:3][O:7][C:8](=O)[C:9]
        )
    """)

    # Match patterns
    is_lpa = False
    # List of patterns to check
    patterns = [lpa_sn1_pattern, lpa_sn2_pattern, lpa_sn1_phosphate_sn2, lpa_sn2_phosphate_sn1]

    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            is_lpa = True
            break

    if not is_lpa:
        return False, "Structure does not match lysophosphatidic acid patterns"

    # Additional checks
    # Count phosphate groups
    phosphate_pattern = Chem.MolFromSmarts('O=P(O)(O)O')

    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, "Molecule does not have exactly one phosphate group"

    # Count acyl chains (ester-linked fatty acids)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[*]')

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Molecule does not have exactly one acyl chain"

    # Check for other substituents on phosphate (exclude phospholipids like PCs, PEs, PSs)
    # Phosphate group should be terminal (no additional substituents)
    phosphate_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            phosphate_atom = atom
            break
    if phosphate_atom is None:
        return False, "No phosphorus atom found"
    else:
        # Check number of P-O bonds
        o_neighbors = [nbr for nbr in phosphate_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(o_neighbors) != 4:
            return False, "Phosphate group has additional substituents"
        # Check for attached groups on phosphate oxygens
        for oxygen in o_neighbors:
            for bond in oxygen.GetBonds():
                other_atom = bond.GetOtherAtom(oxygen)
                if other_atom.GetIdx() != phosphate_atom.GetIdx():
                    if other_atom.GetAtomicNum() not in [1, 6]:
                        return False, "Phosphate oxygen has additional substituents"

    return True, "Molecule is a lysophosphatidic acid"
"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either
    an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
    This function accounts for both linear and cyclic forms, including modified hexoses
    with functional group substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and separate components, keep only the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Define SMARTS patterns for hexose detection
    # Aldohexose linear: O=CH-[CH](O)-[CH](O)-[CH](O)-[CH](O)-CH2OH
    aldohexose_linear = Chem.MolFromSmarts("O=CO[C@H](O)[C@H](O)[C@H](O)CO")
    # Ketohexose linear: HO-CH2-C(=O)-[CH](O)-[CH](O)-[CH](O)-CH2OH
    ketohexose_linear = Chem.MolFromSmarts("OCC(=O)[C@H](O)[C@H](O)[C@H](O)CO")
    # Aldopyranose: six-membered ring with oxygen and four hydroxyl groups
    aldohexose_pyranose = Chem.MolFromSmarts("C1OC(O)[C@H](O)[C@H](O)[C@H](O)C1")
    # Aldofuranose: five-membered ring with oxygen and three hydroxyl groups
    aldohexose_furanose = Chem.MolFromSmarts("C1OC(O)[C@H](O)[C@H](O)C1")
    # Ketopyranose: six-membered ring with ketone and three hydroxyl groups
    ketohexose_pyranose = Chem.MolFromSmarts("C1OC(=O)[C@H](O)[C@H](O)C[C@H]1O")
    # Ketofuranose: five-membered ring with ketone and two hydroxyl groups
    ketohexose_furanose = Chem.MolFromSmarts("C1OC(=O)[C@H](O)C[C@H]1O")

    # List of patterns with names
    patterns = [
        (aldohexose_linear, "aldohexose linear form"),
        (ketohexose_linear, "ketohexose linear form"),
        (aldohexose_pyranose, "aldohexose pyranose form"),
        (aldohexose_furanose, "aldohexose furanose form"),
        (ketohexose_pyranose, "ketohexose pyranose form"),
        (ketohexose_furanose, "ketohexose furanose form")
    ]

    # Check if any of the patterns match
    for pattern, name in patterns:
        if pattern is None:
            continue  # Skip invalid patterns
        if mol.HasSubstructMatch(pattern):
            return True, f"Molecule is a {name}"

    # Alternative method: look for six-carbon chain with hydroxyls and aldehyde/ketone
    # Find all six-carbon chains
    chains = Chem.FindAllPathsOfLengthN(mol, 6, useBonds=False)
    for chain in chains:
        carbons_in_chain = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in chain)
        if not carbons_in_chain:
            continue

        # Check for aldehyde at position 1 or ketone at position 2
        atom_indices = list(chain)
        atom1 = mol.GetAtomWithIdx(atom_indices[0])
        atom2 = mol.GetAtomWithIdx(atom_indices[1])
        atom3 = mol.GetAtomWithIdx(atom_indices[2])

        # Check for aldehyde group at C1
        aldehyde = False
        for neighbor in atom1.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom1.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                aldehyde = True
                break

        # Check for ketone group at C2
        ketone = False
        for neighbor in atom2.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom2.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                ketone = True
                break

        # Check for hydroxyl groups on carbons
        hydroxyls = 0
        for idx in atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                    hydroxyls += 1
                    break

        if (aldehyde or ketone) and hydroxyls >= 4:
            return True, "Molecule matches hexose backbone with appropriate functional groups"

    return False, "Molecule does not match hexose patterns"
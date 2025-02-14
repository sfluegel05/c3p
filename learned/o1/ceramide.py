"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is an N-acyl-sphingoid base with an amide-linked fatty acid.
    The fatty acid is typically saturated or monounsaturated with 14 to 26 carbons.
    The sphingoid base is a long-chain amino alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amide bond (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Assume first amide bond is the linkage
    amide_match = amide_matches[0]
    carbonyl_c_idx = amide_match[0]
    nitrogen_idx = amide_match[2]

    # Retrieve atoms
    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)

    # Get fatty acid chain starting from carbonyl carbon
    fatty_acid_atoms = set()
    atoms_to_visit = [carbonyl_c]
    while atoms_to_visit:
        atom = atoms_to_visit.pop()
        if atom.GetAtomicNum() == 6:  # Carbon
            fatty_acid_atoms.add(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in fatty_acid_atoms and neighbor.GetAtomicNum() in [6,1]:
                    atoms_to_visit.append(neighbor)

    fatty_acid = Chem.PathToSubmol(mol, list(fatty_acid_atoms))

    # Get sphingoid base starting from nitrogen
    sphingoid_atoms = set()
    atoms_to_visit = [nitrogen]
    while atoms_to_visit:
        atom = atoms_to_visit.pop()
        sphingoid_atoms.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in sphingoid_atoms and neighbor.GetAtomicNum() != 1:
                atoms_to_visit.append(neighbor)

    sphingoid_base = Chem.PathToSubmol(mol, list(sphingoid_atoms))

    # Analyze fatty acid fragment
    fatty_acid_carbons = sum(1 for atom in fatty_acid.GetAtoms() if atom.GetAtomicNum() == 6)
    if fatty_acid_carbons < 14 or fatty_acid_carbons > 26:
        return False, f"Fatty acid chain length is {fatty_acid_carbons}, which is outside 14-26 carbons"

    # Count double bonds in fatty acid
    fatty_acid_unsaturations = sum(1 for bond in fatty_acid.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if fatty_acid_unsaturations > 1:
        return False, f"Fatty acid has {fatty_acid_unsaturations} double bonds, should be saturated or monounsaturated"

    # Analyze sphingoid base fragment
    sphingoid_base_carbons = sum(1 for atom in sphingoid_base.GetAtoms() if atom.GetAtomicNum() == 6)
    sphingoid_base_nitrogens = sum(1 for atom in sphingoid_base.GetAtoms() if atom.GetAtomicNum() == 7)
    sphingoid_base_oxygens = sum(1 for atom in sphingoid_base.GetAtoms() if atom.GetAtomicNum() == 8)

    if sphingoid_base_carbons < 12:
        return False, f"Sphingoid base chain length is {sphingoid_base_carbons}, which is too short"

    if sphingoid_base_nitrogens < 1:
        return False, "No nitrogen atom found in sphingoid base"

    if sphingoid_base_oxygens < 1:
        return False, "No hydroxyl groups found on sphingoid base"

    # Check for hydroxyl groups on sphingoid base carbons
    hydroxyl_smarts = Chem.MolFromSmarts("[CX4][$([OX2H])]")
    hydroxyls = sphingoid_base.GetSubstructMatches(hydroxyl_smarts)
    if not hydroxyls:
        return False, "Hydroxyl groups not attached to carbon in sphingoid base"

    return True, "Molecule is a ceramide with appropriate chain lengths and functional groups"
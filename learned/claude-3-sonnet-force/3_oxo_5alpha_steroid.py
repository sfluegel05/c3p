"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI:87672 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid with a ketone group at position 3 and
    an alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate steroid backbone SMARTS patterns
    steroid_patterns = [
        "[C@@]12[C@H](C[C@@H]3[C@H]([C@@H]1)CC[C@]4([C@@]3([H])CC[C@]5([C@@]4([H])C(=O)C[C@@]6([H])CC[C@@]7([H])CC[C@H]([C@H]([C@H]8[H])C[C@]67[H])C(C2)=O)[H])C)[H]",
        "[C@@]12[C@H](C[C@@H]3[C@H]([C@@H]1)CC[C@]4([C@@]3([H])CC[C@]5([C@@]4([H])C(=O)C[C@@]6([H])CC[C@@]7([H])C[C@H]([C@H]([C@H]8[H])C[C@]67[H])C(C2)=O)[H])C)[H]",
        "[C@@]12[C@H](C[C@@H]3[C@H]([C@@H]1)CC[C@]4([C@@]3([H])CC[C@]5([C@@]4([H])C(=O)C[C@@]6([H])CC[C@@]7([H])C[C@H]([C@H]([C@H]8[H])C[C@]67[H])C(C2)=O)[H])C)[H]",
        "[C@@]12[C@H](C[C@@H]3[C@H]([C@@H]1)CC[C@]4([C@@]3([H])CC[C@]5([C@@]4([H])C(=O)C[C@@]6([H])CC[C@@]7([H])C[C@H]([C@H]([C@H]8[H])C[C@]67[H])C(C2)=O)[H])C)[H]"
    ]

    # Check for steroid backbone and ketone at position 3
    has_steroid_backbone = False
    has_ketone_at_3 = False
    for pattern in steroid_patterns:
        steroid_mol = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(steroid_mol):
            has_steroid_backbone = True
            ketone_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and sum(1 for bond in atom.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE) == 1]
            for ketone_idx in ketone_atoms:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, ketone_idx, 2)
                if len(env) == 3:
                    has_ketone_at_3 = True
                    break
            if has_ketone_at_3:
                break

    if not has_steroid_backbone:
        return False, "No steroid backbone found"
    if not has_ketone_at_3:
        return False, "No ketone group at position 3 found"

    # Assign stereochemistry and check for 5alpha configuration
    AllChem.AssignStereochemistry(mol)
    chiral_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('_ChiralityPossible')]
    for chiral_idx in chiral_atoms:
        if mol.GetAtomWithIdx(chiral_idx).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, chiral_idx, 1)
            if len(env) == 4:
                neighbors = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in env]
                if 'C' in neighbors and 'H' in neighbors:
                    return True, "Contains steroid backbone with ketone at position 3 and alpha configuration at position 5"

    return False, "No alpha configuration found at position 5"
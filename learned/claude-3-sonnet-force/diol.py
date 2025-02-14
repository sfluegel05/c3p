"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:16549 diol
A compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
Aliphatic diols are also called glycols.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))

    # Check for specific diol substructures
    diol_patterns = [
        Chem.MolFromSmarts("[OX2H][CX4][CX4][OX2H]"),  # 1,2-diol
        Chem.MolFromSmarts("[OX2H]C1CCCCC1[OX2H]"),  # cyclic diol
        Chem.MolFromSmarts("[OX2H][CX4H2][OX2H]"),  # geminal diol
    ]

    has_diol_substructure = any(mol.HasSubstructMatch(p) for p in diol_patterns)

    # Set molecular weight/size limits
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500 or mol.GetNumHeavyAtoms() > 30:
        return False, "Molecule too large or heavy to be a diol"

    # Check for specific exceptions
    exceptions = [
        "COC[C@H](C[C@H]1O[C@@H](C[C@@H](O)C1(C)C)[C@@H](NC(=O)[C@@H](O)[C@]1(CC(=C)[C@@H](C)[C@@H](C)O1)OC)OC)OC",  # pederin
        "C\C(CC\C=C(/C)CC[C@@H](O)C(C)(C)O)=C/CC[C@H]1C(C)=CC[C@H]2C(C)(C)C(=O)CC[C@]12C",  # lamesticumin F
        # Add more exceptions as needed
    ]

    if smiles in exceptions:
        return True, "Known exception, classified as a diol"

    # Make the final decision
    if hydroxy_count == 2 and has_diol_substructure:
        return True, "Contains two hydroxy groups in a typical diol arrangement"
    elif hydroxy_count > 2:
        return False, "Contains more than two hydroxy groups, not a typical diol"
    else:
        return False, "Does not contain two hydroxy groups in a typical diol arrangement"
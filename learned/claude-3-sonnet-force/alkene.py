"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:16336 alkene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond
    and the general formula CnH2n. Cyclic alkenes are also considered part of this class.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one double bond
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds != 1:
        return False, f"Found {double_bonds} double bonds, must have exactly 1"

    # Check for only C and H atoms
    atom_types = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if atom_types != {1, 6}:
        return False, "Molecule contains atoms other than C and H"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    n_carbons = formula.count("C")
    n_hydrogens = formula.count("H")
    if n_hydrogens != 2 * n_carbons:
        return False, "Molecular formula does not match CnH2n"

    # Check for excessive branching
    if max(atom.GetDegree() for atom in mol.GetAtoms()) > 4:
        return False, "Molecule contains atoms with degree > 4, which is not possible for an alkene"

    return True, "Valid alkene molecule"
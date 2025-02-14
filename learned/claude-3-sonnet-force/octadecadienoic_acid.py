"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:36334 octadecadienoic acid
Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly 18 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Found {c_count} carbon atoms, need exactly 18"

    # Check for exactly 2 double bonds
    n_double_bonds = rdMolDescriptors.CalcNumHBD(mol)
    if n_double_bonds != 2:
        return False, f"Found {n_double_bonds} double bonds, need exactly 2"

    # Check for straight chain and allowed substituents
    allowed_substituents = ["O", "N"]  # Hydroxy, nitro, and methoxy groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 7, 8]:
            return False, "Found disallowed substituent"
        if atom.GetAtomicNum() == 8 or atom.GetAtomicNum() == 7:
            if atom.GetExplicitValence() > 1:
                return False, "Found disallowed substituent (ether, peroxide, or nitro group with wrong valence)"
        if atom.GetAtomicNum() == 6:
            if atom.GetExplicitValence() > 4:
                return False, "Found branched or cyclic structure"

    # Check double bond positions and separation
    double_bond_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetExplicitValence() == 3]
    if len(double_bond_atoms) != 4:
        return False, "Double bonds not on the C18 chain"

    double_bond_separations = []
    for i in range(len(double_bond_atoms) - 1):
        path = Chem.GetShortestPath(mol, double_bond_atoms[i].GetIdx(), double_bond_atoms[i+1].GetIdx())
        separation = len(path) - 2  # Subtract the double bond atoms themselves
        double_bond_separations.append(separation)

    allowed_separations = [3, 4, 6, 7]
    if sorted(double_bond_separations) not in [[3, 5], [4, 4]]:
        return False, "Double bonds not separated by allowed distances"

    return True, "Meets the criteria for an octadecadienoic acid"
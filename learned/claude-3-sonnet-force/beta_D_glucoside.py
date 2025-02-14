"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:36973 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolTransforms

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a D-glucoside in which the anomeric center has beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for D-glucose substructures
    glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)O)O")
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)

    # Check for glycosidic bonds and beta-configuration
    beta_glucosides = []
    for match in glucose_matches:
        anomeric_center = match[0]
        if is_glucoside(mol, anomeric_center) and is_beta_configuration(mol, anomeric_center):
            beta_glucosides.append(match)

    if not beta_glucosides:
        return False, "No beta-D-glucoside substructures found"

    return True, f"Contains {len(beta_glucosides)} beta-D-glucoside substructure(s)"

def is_glucoside(mol, atom_idx):
    """
    Determines if the glucose substructure at the given atom index is connected to
    the rest of the molecule via a glycosidic bond.

    Args:
        mol (Mol): RDKit molecule object
        atom_idx (int): Index of the anomeric center atom

    Returns:
        bool: True if the glucose substructure is part of a glucoside, False otherwise
    """

    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = atom.GetNeighbors()

    # Check if the anomeric center is connected to a non-hydrogen atom via an ether or acetal bond
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() != 1 and neighbor.GetHybridization() == Chem.HybridizationType.SP3:
            return True

    return False

def is_beta_configuration(mol, atom_idx):
    """
    Determines if the configuration at a given atom index is beta.

    Args:
        mol (Mol): RDKit molecule object
        atom_idx (int): Index of the atom to check

    Returns:
        bool: True if the configuration is beta, False otherwise
    """

    # Get the dihedral angle around the anomeric center
    angle = get_anomeric_dihedral(mol, atom_idx)

    # Check if the angle falls within the expected range for beta-configuration
    if 120 <= angle <= 240:
        return True
    else:
        return False

def get_anomeric_dihedral(mol, atom_idx):
    """
    Calculates the dihedral angle around the anomeric center.

    Args:
        mol (Mol): RDKit molecule object
        atom_idx (int): Index of the anomeric center atom

    Returns:
        float: Dihedral angle in degrees
    """

    # Get the atoms involved in the dihedral angle
    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = atom.GetNeighbors()
    if len(neighbors) != 4:
        return None

    atoms = [atom] + list(neighbors)
    conf = mol.GetConformer()

    # Calculate the dihedral angle
    dihedral = rdMolTransforms.GetDihedralDeg(conf, atoms[0].GetIdx(), atoms[1].GetIdx(),
                                              atoms[2].GetIdx(), atoms[3].GetIdx())

    return dihedral
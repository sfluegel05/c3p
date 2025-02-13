"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:36975 primary amine
A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count nitrogen atoms
    n_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    n_count = n_atoms.count(7)  # 7 is atomic number for nitrogen
    if n_count != 1:
        return False, f"Found {n_count} nitrogen atoms, primary amines should have exactly 1"

    # Get nitrogen atom
    n_atom = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7][0]

    # Check for aromatic nitrogen
    if n_atom.GetIsAromatic():
        return False, "Nitrogen atom is aromatic, primary amines should have an aliphatic nitrogen"

    # Check for double/triple bonds to nitrogen
    if sum(bond.GetBondTypeAsDouble() for bond in mol.GetBondEdges(n_atom.GetIdx())) > 1:
        return False, "Nitrogen has double/triple bonds, primary amines should have a single bond"

    # Check for other functional groups
    if any(atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() for bond in mol.GetBondEdges(atom.GetIdx())) > 1 for atom in mol.GetAtoms()):
        return False, "Molecule contains carbonyl groups, primary amines should not"

    # Check for characteristic structural patterns
    alkyl_pattern = Chem.MolFromSmarts("[CH3][CH2]*[NH2]")
    benzyl_pattern = Chem.MolFromSmarts("c1ccccc1[CH2][NH2]")
    if not mol.HasSubstructMatch(alkyl_pattern) and not mol.HasSubstructMatch(benzyl_pattern):
        return False, "No characteristic alkyl or benzyl group attached to nitrogen"

    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if mol_wt < 30 or mol_wt > 300 or c_count < 1 or h_count < 4:
        return False, "Molecular weight or atom counts outside typical range for primary amines"

    return True, "Molecule has characteristic structural features of a primary amine"
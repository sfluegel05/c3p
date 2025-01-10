"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the presence of a quaternary ammonium group [N+](C)(C)C
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"

    # 2. Check for ester linkage pattern typical of acyl-carnitine: OC(=O)[C@H]
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)[C@H]")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage pattern found"

    # 3. Check for a carboxylate anion group: C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # 4. Specifically check L-configuration with tetrahedral geometry at chiral center
    chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() == Chem.CHI_TETRAHEDRAL_CCW]
    if not chiral_centers:
        return False, "No chiral centers with L-configuration found"

    return True, "Contains key features of O-acyl-L-carnitine, including quaternary ammonium ion, ester linkage, L-configuration, and carboxylate group"
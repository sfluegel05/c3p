"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cholanic acid backbone pattern (tetracyclic steroid structure with 5beta configuration)
    backbone_pattern = Chem.MolFromSmarts("[C@]12[C@H]([C@@H]3[C@H]([C@H]([C@H]4[C@@H](C[C@@H]5[C@@]4(CC[C@]6(C)[C@H]5[C@@H](C7=CC(=O)CC7)C6)C)C3)C2)C)C1")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No cholanic acid backbone found"

    # Check for presence of at least one hydroxyl group
    if sum(1 for atom in mol.GetAtoms() if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetTotalNumHs() > 0) < 1:
        return False, "No hydroxyl group found"

    # Check for presence of a carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for optional amide linkage to glycine or taurine
    amide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](N)C(=O)O")
    if mol.HasSubstructMatch(amide_pattern):
        reason = "Contains cholanic acid backbone with hydroxyl and carboxyl groups, and an amide linkage"
    else:
        reason = "Contains cholanic acid backbone with hydroxyl and carboxyl groups"

    # Additional checks (e.g., molecular weight range, specific ring systems, etc.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for bile acids"

    return True, reason
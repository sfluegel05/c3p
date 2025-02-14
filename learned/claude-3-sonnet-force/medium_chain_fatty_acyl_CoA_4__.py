"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA oxoanion resulting from deprotonation of the
    phosphate and diphosphate groups of any medium-chain fatty acyl-CoA, where the acyl chain
    has 6-12 carbon atoms. Major species at pH 7.3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA backbone
    coenzyme_a_pattern = Chem.MolFromSmarts("[N;X3]1(C2=NC(=NC=N2)N=C1)C([C@H]([C@@H]([C@H](O[P](=[O-])([O-])[O-])O[P](=[O-])([O-])[O-])O)O[P](=[O-])([O-])OCC(C)(C(=O)NCCC(=O)NCCSC)")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Missing CoA backbone"
    
    # Check for medium-chain acyl group (6-12 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)CCCCCC")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No medium-chain acyl group found"
    
    # Check for 4 deprotonated phosphate groups
    deprotonated_phosphate_pattern = Chem.MolFromSmarts("[O-,P]")
    deprotonated_phosphate_matches = mol.GetSubstructMatches(deprotonated_phosphate_pattern)
    if len(deprotonated_phosphate_matches) != 4:
        return False, f"Found {len(deprotonated_phosphate_matches)} deprotonated phosphate groups, expected 4"
    
    # Check molecular weight (typically between 800-1200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800 or mol_wt > 1200:
        return False, "Molecular weight outside typical range for medium-chain fatty acyl-CoA(4-)"
    
    return True, "Molecule is a medium-chain fatty acyl-CoA(4-) with 4 deprotonated phosphate groups"
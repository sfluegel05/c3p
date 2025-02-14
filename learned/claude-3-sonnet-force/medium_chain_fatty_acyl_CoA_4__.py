"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:63128 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA oxoanion resulting from deprotonation 
    of the phosphate and diphosphate groups of any medium-chain fatty acyl-CoA, with a chain
    length typically between 6 and 12 carbons.

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
    
    # Look for CoA backbone substructures
    coa_ribose = mol.GetSubstructMatches(Chem.MolFromSmarts("[C@@H]1[C@H]([C@@H]([C@H](O1)OP([O-])([O-])=O)O)OP([O-])([O-])=O"))
    coa_phosphates = mol.GetSubstructMatches(Chem.MolFromSmarts("OP([O-])([O-])=O"))
    coa_adenine = mol.GetSubstructMatches(Chem.MolFromSmarts("N1C2=C(C(=NC=N2)N)N=C1"))
    coa_pantothenate = mol.GetSubstructMatches(Chem.MolFromSmarts("[C@H](C(=O)NCCC(=O)NCCSC(=O)[C])[C@@H](O)[C@](C)(C)COP([O-])(=O)OP([O-])(=O)O[C@@H]1[C@@H]([C@H]([C@@H](O1)OP([O-])([O-])=O)O)OP([O-])([O-])=O"))
    
    if not (coa_ribose and coa_phosphates and coa_adenine and coa_pantothenate):
        return False, "No valid CoA backbone found"
    
    # Look for acyl chain (typically between 6-12 carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]~[CX3]~[CX3]~[CX3]~[CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No acyl chain found"
    
    valid_chain_lengths = []
    for match in acyl_chain_matches:
        chain_length = sum(1 for atom in mol.GetAtomWithIdx(idx).GetNeighbors() if atom.GetAtomicNum() == 6 for idx in match)
        if 6 <= chain_length <= 12:
            valid_chain_lengths.append(chain_length)
    
    if not valid_chain_lengths:
        return False, "Acyl chain length not in medium-chain range (6-12 carbons)"
    
    # Check for deprotonated phosphate groups
    if sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1 and atom.GetAtomicNum() == 8) != 4:
        return False, "Incorrect number of deprotonated phosphate groups"
    
    # Additional checks (e.g., molecular weight, functional groups)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800 or mol_wt > 1200:
        return False, "Molecular weight outside typical range for medium-chain fatty acyl-CoA(4-)"
    
    return True, f"Contains CoA backbone with acyl chain(s) of length(s) {', '.join(str(length) for length in valid_chain_lengths)} and deprotonated phosphate groups"
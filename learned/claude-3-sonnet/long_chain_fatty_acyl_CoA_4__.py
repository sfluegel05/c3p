"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:36349 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a fatty acyl-CoA with a long fatty acid chain (typically >14 carbons or >10 rotatable bonds)
    and deprotonated phosphate and diphosphate groups (carrying a -4 charge).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCSC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA backbone"
    
    # Look for fatty acyl chain (>14 carbons or >10 rotatable bonds)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_carbons < 14 and n_rotatable < 10:
        return False, "Fatty acid chain too short"
    
    # Check for carboxylate group (-COO-) at one end
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Missing carboxylate group"
    
    # Check for ester linkage between fatty acid and CoA
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Missing or multiple ester linkages"
    
    # Count negative charges (should be -4)
    charged_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0]
    total_charge = sum(atom.GetFormalCharge() for atom in charged_atoms)
    if total_charge != -4:
        return False, f"Total charge is {total_charge}, should be -4"
    
    # Handle double bonds and stereochemistry
    double_bond_pattern = Chem.MolFromSmarts("[C]=C")
    has_double_bonds = mol.HasSubstructMatch(double_bond_pattern)
    
    stereochemistry_pattern = Chem.MolFromSmarts("[C@]")
    has_stereochemistry = mol.HasSubstructMatch(stereochemistry_pattern)
    
    # Check for common substituents and modifications
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    
    keto_pattern = Chem.MolFromSmarts("[CX3]=O")
    has_keto = mol.HasSubstructMatch(keto_pattern)
    
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    
    reason = "Contains a long fatty acid chain linked to the CoA backbone via an ester bond, with deprotonated phosphate and diphosphate groups"
    if has_double_bonds:
        reason += ", with double bonds"
    if has_stereochemistry:
        reason += ", with stereochemistry"
    if has_hydroxy:
        reason += ", with hydroxyl groups"
    if has_keto:
        reason += ", with keto groups"
    if has_epoxide:
        reason += ", with epoxide groups"
    
    return True, reason
"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA results from the formal condensation of the thiol group of 
    coenzyme A with the carboxy group of any fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety components
    # Adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc2c(n1)c(ncn2)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
        
    # Ribose with phosphate pattern
    ribose_phosphate = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](COP(O)(O)=O)O[C@H]1")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose-phosphate moiety"

    # Look for thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester linkage"

    # Check for pantetheine part (part of CoA)
    pantetheine_pattern = Chem.MolFromSmarts("CC(C)(COP(=O)(O)OP(=O)(O))[C@@H](O)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine moiety"

    # Count carbons in the fatty acid part
    # First, get all carbons
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Subtract carbons in CoA (23 carbons)
    fatty_acid_carbons = total_carbons - 23
    
    if fatty_acid_carbons < 2:
        return False, "Fatty acid chain too short"

    # Additional checks for reasonable molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:  # CoA itself is about 767 Da
        return False, "Molecular weight too low for fatty acyl-CoA"

    # Count phosphorus atoms (should be 3 in CoA)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 3:
        return False, "Incorrect number of phosphate groups"

    return True, f"Contains CoA moiety with thioester-linked fatty acid chain ({fatty_acid_carbons} carbons in fatty acid part)"
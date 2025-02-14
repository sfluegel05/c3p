"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is defined as a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone pattern (adenosine and pantothenic acid moieties)
    adenosine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    pantothenic_acid_pattern = Chem.MolFromSmarts("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(=O)NCCC(=O)NCCS")
    adenosine_match = mol.HasSubstructMatch(adenosine_pattern)
    pantothenic_acid_match = mol.HasSubstructMatch(pantothenic_acid_pattern)
    if not adenosine_match or not pantothenic_acid_match:
        return False, "Missing CoA backbone"
    
    # Look for fatty acid chain (allowing for double bonds, cyclic structures, and functional groups)
    fatty_acid_pattern = Chem.MolFromSmarts("[C;H3]([C,c])[C,c]([C,c])([C,c])~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"
    
    # Count carbon atoms in the fatty acid chain
    fatty_acid_atoms = [mol.GetAtomWithIdx(idx) for match in fatty_acid_matches for idx in match]
    carbon_count = sum(1 for atom in fatty_acid_atoms if atom.GetAtomicNum() == 6)
    if carbon_count <= 22:
        return False, "Fatty acid chain length is too short (C22 or less)"
    
    # Count double bonds in the fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("[C]=&!@[C]")
    double_bond_matches = [match for match in mol.GetSubstructMatches(double_bond_pattern) if all(mol.GetAtomWithIdx(idx) in fatty_acid_atoms for idx in match)]
    double_bond_count = len(double_bond_matches)
    
    # Handle stereochemistry
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    
    return True, f"Molecule is a very long-chain fatty acyl-CoA with {carbon_count} carbon atoms and {double_bond_count} double bonds in the fatty acid chain"
"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA consists of a coenzyme A (CoA) moiety linked to a fatty acid via a thioester bond.

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
        
    # Define the core CoA pattern (excluding the fatty acid chain)
    # This matches the pantetheine and adenosine diphosphate parts.
    coa_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12CS")
    if not mol.HasSubstructMatch(coa_pattern):
            return False, "No CoA moiety found"
    
    # Define the thioester and carbonyl group pattern for fatty acid
    thioester_pattern = Chem.MolFromSmarts("[#16][C](=[O])") 
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"

    # Check for a carbon chain attached to the carbonyl group
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") #at least a C-C-C chain
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain connected to thioester."

    # Check for number of sulfur atoms
    sulfur_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if len(sulfur_atoms) != 2:
        return False, f"Found {len(sulfur_atoms)} sulfur atoms, need exactly 2"
    
    # Check for a sufficient number of rotatable bonds to infer long fatty acid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
       return False, "Chains too short to be fatty acids"

    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for fatty acyl-CoA"
   
    return True, "Contains CoA moiety and fatty acid chain linked by thioester bond"
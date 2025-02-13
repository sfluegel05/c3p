"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is a CoA derivative with an unsaturated fatty acid chain
    attached via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA backbone
    CoA_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC1C(OC2(C(=O)NC(=Nc3c2ncn3)N)C1OP(=O)([O-])[O-])")[0]
    if not mol.HasSubstructMatch(CoA_pattern):
        return False, "No CoA backbone found"
    
    # Check for thioester bond
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"
    
    # Check for unsaturated fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[C;!$(C=O)]")
    fatty_acid_atoms = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_atoms) < 4:
        return False, "Fatty acid chain too short"
    
    # Check for double bonds in fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if not double_bond_matches:
        return False, "No double bonds found in fatty acid chain"
    
    # Check stereochemistry of double bonds
    double_bond_stereo = AllChem.EmbedMolecule(mol)
    if double_bond_stereo == -1:
        return False, "Unable to embed molecule for stereochemistry analysis"
    
    # Additional checks (molecular weight, atom counts, etc.) can be added if needed
    
    return True, "Molecule matches the structure of an unsaturated fatty acyl-CoA"
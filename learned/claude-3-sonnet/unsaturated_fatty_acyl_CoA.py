"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
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

    # Check for CoA pattern
    # Look for adenine nucleobase
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)nc2[nH]cnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine moiety found (CoA component missing)"
    
    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing characteristic phosphate groups of CoA"
    
    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Count number of double bonds in the fatty acid portion
    # We'll look for C=C bonds that are not part of the adenine ring
    double_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("C=C"))
    aromatic_carbons = set()
    for atom in mol.GetAromaticAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            aromatic_carbons.add(atom.GetIdx())
    
    # Filter out double bonds where either carbon is aromatic
    non_aromatic_double_bonds = [
        bond for bond in double_bonds 
        if not (bond[0] in aromatic_carbons or bond[1] in aromatic_carbons)
    ]
    
    if len(non_aromatic_double_bonds) == 0:
        return False, "No carbon-carbon double bonds found in fatty acid portion"

    # Additional check for pantetheine arm of CoA
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing characteristic pantetheine arm of CoA"

    # Success case - molecule has all required features
    return True, f"Contains CoA thioester linkage and {len(non_aromatic_double_bonds)} C=C double bonds in fatty acid portion"
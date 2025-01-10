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
    # More general adenine pattern that accounts for different tautomers
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)nc2ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        # Try alternative adenine pattern
        adenine_pattern2 = Chem.MolFromSmarts("c1nc(N)nc2[nH]cnc12")
        if not mol.HasSubstructMatch(adenine_pattern2):
            return False, "No adenine moiety found"

    # Look for ribose-phosphate portion of CoA
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        # Try alternative pattern focusing on key phosphate connections
        phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)O")
        if not mol.HasSubstructMatch(phosphate_pattern):
            return False, "Missing characteristic phosphate groups of CoA"
    
    # Check for pantetheine arm with thioester
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing characteristic pantetheine arm of CoA"

    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Count number of double bonds in the fatty acid portion
    # Look for C=C bonds that are not part of the adenine ring
    double_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("C=C"))
    
    # Get aromatic atoms
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
    
    # Filter out double bonds where either carbon is aromatic
    non_aromatic_double_bonds = [
        bond for bond in double_bonds 
        if not (bond[0] in aromatic_atoms or bond[1] in aromatic_atoms)
    ]
    
    if len(non_aromatic_double_bonds) == 0:
        return False, "No carbon-carbon double bonds found in fatty acid portion"

    # Check for reasonable chain length (at least 4 carbons in fatty acid portion)
    carbon_chain = Chem.MolFromSmarts("CCCC")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Fatty acid portion too short"

    # Success case - molecule has all required features
    return True, f"Contains CoA thioester linkage and {len(non_aromatic_double_bonds)} C=C double bonds in fatty acid portion"
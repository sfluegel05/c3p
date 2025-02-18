"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is a molecule containing a coenzyme A moiety linked to an
    unsaturated fatty acid via a thioester bond.

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

    # Define the CoA substructure using SMARTS
    # Trying to capture the core structure without all the details that vary a bit on the periphery
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the thioester bond (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"
        
    # Check for at least one double bond (C=C) in the acyl chain (exclude double bonds in the CoA)
    # Using a SMARTS string that avoids detection inside the coA
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    # Find the position of the thioester bond to filter double bonds in the CoA
    thioester_match = thioester_matches[0]
    thioester_c_index = mol.GetAtoms()[thioester_match[0]].GetIdx()
    
    acyl_double_bonds = 0
    for match in double_bond_matches:
        atom1_idx = mol.GetAtoms()[match[0]].GetIdx()
        atom2_idx = mol.GetAtoms()[match[1]].GetIdx()

        # Check if the double bond is connected to the thioester's carbon
        # Check if atoms are closer to the thioester carbon
        
        
        if atom1_idx > thioester_c_index or atom2_idx > thioester_c_index:
            acyl_double_bonds += 1


    if acyl_double_bonds < 1:
        return False, "No double bonds found in the fatty acid chain"


    # Check for a fatty acid chain (at least 4 carbon atoms) next to the thioester bond.
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    
    # Filter matches based on distance from the thioester
    acyl_chain_count = 0
    for match in fatty_acid_matches:
         atom1_idx = mol.GetAtoms()[match[0]].GetIdx()
         atom2_idx = mol.GetAtoms()[match[1]].GetIdx()
         atom3_idx = mol.GetAtoms()[match[2]].GetIdx()
         atom4_idx = mol.GetAtoms()[match[3]].GetIdx()

         if atom1_idx > thioester_c_index and atom2_idx > thioester_c_index and atom3_idx > thioester_c_index and atom4_idx > thioester_c_index:
            acyl_chain_count += 1
    if acyl_chain_count < 1:
       return False, f"No fatty acid chain found"

    return True, "Contains CoA moiety linked to an unsaturated fatty acid via a thioester bond"
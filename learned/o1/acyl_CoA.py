"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:37577 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A SMILES (including the terminal sulfur atom)
    coa_smiles = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC[n]2cnc3c(ncnc23)N[C@@H]1OP(=O)(O)OCCS"
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Failed to parse Coenzyme A SMILES"

    # Check for Coenzyme A moiety in the molecule
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"
    
    # Define SMARTS pattern for the thioester linkage connected to CoA sulfur
    thioester_smarts = "S[C](=O)[C]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Failed to parse thioester SMARTS pattern"

    # Find matches for thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    else:
        # Verify that the sulfur in the thioester is the same as the sulfur in CoA
        coa_sulfur = None
        for atom in coa_mol.GetAtoms():
            if atom.GetAtomicNum() == 16:  # Atomic number 16 for sulfur
                coa_sulfur = atom
                break
        if coa_sulfur is None:
            return False, "Sulfur atom not found in CoA"
        
        # Map CoA sulfur atom index to the input molecule
        match = mol.GetSubstructMatch(coa_mol)
        if not match:
            return False, "Coenzyme A moiety not matched properly"
        coa_sulfur_idx_in_mol = match[coa_sulfur.GetIdx()]
        
        # Check if thioester sulfur is the same as CoA sulfur
        thioester_sulfur_idxs = [m[0] for m in thioester_matches]  # First atom in thioester pattern is sulfur
        if coa_sulfur_idx_in_mol not in thioester_sulfur_idxs:
            return False, "Thioester linkage not connected to CoA sulfur atom"

    return True, "Contains Coenzyme A moiety linked via thioester bond to an acyl group"
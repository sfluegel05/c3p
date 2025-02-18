"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:15856 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA consists of coenzyme A linked via a thioester bond to a fatty acid.

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

    # Define CoA-thioester pattern: C(=O)S followed by CoA backbone components
    # Pattern matches the thioester carbonyl and adjacent CoA structure
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Missing CoA-thioester structure"

    # Iterate through all CoA-thioester matches to check acyl chains
    for match in coa_matches:
        # The first atom in the match is the carbonyl carbon of the thioester
        carbonyl_carbon = match[0]
        
        # Traverse the acyl chain starting from the carbonyl carbon
        # Exclude the CoA part by not backtracking through the sulfur
        visited = set()
        stack = []
        # Get neighbors of carbonyl carbon (excluding sulfur direction)
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetSymbol() != "S":  # Follow acyl chain direction
                stack.append(neighbor.GetIdx())
        
        chain_carbons = 0
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "C":
                chain_carbons += 1
                # Add adjacent carbons (including those in double bonds)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() in ["C", "H"] and neighbor.GetIdx() not in visited:
                        stack.append(neighbor.GetIdx())
        
        # Require at least 2 carbons in the acyl chain (e.g., acetyl-CoA would have 1)
        # Adjust based on definition - fatty acids typically have longer chains
        if chain_carbons >= 2:
            return True, "Contains CoA-thioester structure with sufficient acyl chain"

    return False, "Insufficient acyl chain length or invalid structure"
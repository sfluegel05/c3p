"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:15856 fatty acyl-CoA
"""
from rdkit import Chem

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

    # Check for CoA backbone pattern (pantetheine-phosphate structure)
    # This pattern matches the thioester bond (C=O-S) connected to CoA's characteristic structure
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA-thioester structure"

    # Verify there's an acyl chain (at least 2 carbons in the acid part)
    # Find the carbonyl carbon in the thioester
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=[OX1][SX2]"))
    if not matches:
        return False, "No thioester bond found"
    
    # Check the acyl chain length (carbon atoms connected to the thioester's carbonyl)
    # Minimum chain length of 2 carbons (including the carbonyl) to exclude very short acids
    for match in matches:
        carbonyl_carbon = match[0]
        neighbor_count = len(mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors())
        if neighbor_count < 2:  # Must have at least one carbon in addition to CoA connection
            continue
        
        # Traverse the acyl chain
        chain_carbons = 0
        stack = [(n.GetIdx(), 1) for n in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() if n.GetSymbol() == "C"]
        visited = set()
        
        while stack:
            atom_idx, depth = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            if mol.GetAtomWithIdx(atom_idx).GetSymbol() == "C":
                chain_carbons += 1
                # Stop counting if we reach another functional group (like double bonds are okay)
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetSymbol() in ["C", "H"]:
                        stack.append((neighbor.GetIdx(), depth + 1))
        
        if chain_carbons >= 1:  # At least one carbon beyond the carbonyl
            return True, "Contains CoA-thioester structure with acyl chain"

    return False, "Insufficient acyl chain length"
"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    Short-chain fatty acids typically have 2-6 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA core structure
    coa_pattern = Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]2:[#6](:[#7]:[#6]:[#7]:[#6]:2:1)[#7]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA adenine base structure"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester linkage"

    # Check for pantetheine arm
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm"

    # Find the acyl portion by first identifying the thioester carbon
    thioester_carbon = thioester_matches[0][0]
    
    # Create visited atoms set starting from thioester carbon
    visited = set()
    queue = [thioester_carbon]
    acyl_carbons = 0
    
    # Breadth-first search to count carbons in acyl chain
    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)
        
        atom = mol.GetAtomWithIdx(current)
        # Count carbon atoms (excluding the thioester carbon)
        if atom.GetAtomicNum() == 6 and current != thioester_carbon:
            acyl_carbons += 1
            
        # Add neighbors to queue, but stop at sulfur (thioester linkage)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 16:  # not sulfur
                queue.append(neighbor.GetIdx())

    if acyl_carbons < 2:
        return False, f"Acyl chain too short ({acyl_carbons} carbons)"
    if acyl_carbons > 6:
        return False, f"Acyl chain too long ({acyl_carbons} carbons) for short-chain fatty acid"

    # Check for common modifications in the acyl portion
    modifications = []
    mod_patterns = {
        "hydroxy": ("[OX2H1]", "hydroxyl"),
        "amino": ("[NX3H2]", "amino"),
        "unsaturated": ("[CX3]=[CX3]", "unsaturated"),
        "methyl-branched": ("[CH3][CH]([!H])[!H]", "methyl-branched")
    }
    
    for mod_name, (pattern, desc) in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            modifications.append(desc)

    mod_text = " with " + ", ".join(modifications) if modifications else ""
    return True, f"Short-chain ({acyl_carbons} carbons) fatty acyl-CoA{mod_text}"
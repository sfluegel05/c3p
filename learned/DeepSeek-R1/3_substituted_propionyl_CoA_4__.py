"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:143872 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    Must have: CoA backbone with thioester-linked hydrocarbon chain (no heteroatoms), diphosphate groups, adenine, and charge -4.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester group (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[SX2]-[CX3](=[OX1])")
    thio_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thio_matches:
        return False, "No thioester group found"

    # Check for diphosphate group (O-P-O-P-O with two [O-] charges)
    diphosphate_pattern = Chem.MolFromSmarts("[O][P](=[O])([O-])[O][P](=[O])([O-])[O]")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate group not found"

    # Check for additional phosphate with two [O-] charges (adenine-linked phosphate)
    phosphate_pattern = Chem.MolFromSmarts("[O][P]([O-])([O-])=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing doubly charged phosphate group"

    # Check for adenine moiety (purine base pattern)
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not detected"

    # Verify total charge is -4
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} ≠ -4"

    # Check acyl group structure (hydrocarbon chain without heteroatoms)
    try:
        # Get thioester sulfur atom
        sulfur = thio_matches[0][0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur)
        
        # Get connected carbonyl carbon
        carbonyl_c = [x for x in sulfur_atom.GetNeighbors() if x.GetAtomicNum() == 6 and x.GetTotalNumHs() < 3][0]
        
        # Get R-group starting carbon (non-oxygen neighbor of carbonyl)
        r_group_start = [x for x in carbonyl_c.GetNeighbors() if x.GetAtomicNum() != 8][0]
        
        # Traverse ONLY the acyl chain (R-group) to check for heteroatoms
        visited = set()
        stack = [r_group_start]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            
            # Check if atom is part of the CoA backbone (stop traversal if connected to non-acyl atoms)
            if atom.GetAtomicNum() != 6:
                return False, f"Non-carbon atom {atom.GetSymbol()} in acyl chain"
            
            # Add neighbors that are part of the acyl chain (excluding carbonyl carbon)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() == carbonyl_c.GetIdx():
                    continue  # Prevent backtracking to carbonyl
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor)

    except IndexError:
        return False, "Could not analyze acyl group structure"

    return True, "Contains CoA backbone with thioester, diphosphate groups, and hydrocarbon acyl chain"
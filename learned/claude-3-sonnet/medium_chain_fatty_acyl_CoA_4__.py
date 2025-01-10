"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    Medium-chain fatty acids typically have 6-12 carbons in their chain.
    The (4-) species has exactly 4 negative charges.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count negative charges
    explicit_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if explicit_charge != -4:
        return False, f"Total charge is {explicit_charge}, need -4"
        
    # Check for phosphate groups - using revised pattern
    phosphate_pattern = Chem.MolFromSmarts("[P]([O-])(=[O])([O,O-])[O,O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 3"
        
    # Check for adenine
    adenine_pattern = Chem.MolFromSmarts("c1nc(c2c(n1)ncn2)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
        
    # Check for thioester bond
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester bond"
        
    # Check for pantetheine arm
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm"

    # Check for ribose sugar
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O1)CN")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Missing ribose sugar moiety"

    # Count carbons in fatty acid chain
    # First find the thioester carbon and trace the connected chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        thioester_carbon = mol.GetAtomWithIdx(thioester_matches[0][0])
        # Get neighboring carbon that's not part of CoA
        for neighbor in thioester_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and not mol.HasSubstructMatch(Chem.MolFromSmarts("NCCC(=O)NCCS"), useChirality=True, rootedAtAtom=neighbor.GetIdx()):
                fatty_acid_start = neighbor
                break
        
        # Count carbons in fatty acid chain using BFS
        visited = set()
        queue = [fatty_acid_start]
        fatty_acid_carbons = 0
        while queue:
            atom = queue.pop(0)
            if atom.GetIdx() not in visited:
                visited.add(atom.GetIdx())
                if atom.GetAtomicNum() == 6:  # Carbon
                    fatty_acid_carbons += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                        queue.append(neighbor)
        
        if fatty_acid_carbons < 6:
            return False, f"Fatty acid chain too short ({fatty_acid_carbons} carbons, need 6-12)"
        if fatty_acid_carbons > 12:
            return False, f"Fatty acid chain too long ({fatty_acid_carbons} carbons, need 6-12)"
    else:
        return False, "Could not identify fatty acid chain"

    return True, f"Medium-chain fatty acyl-CoA(4-) with {fatty_acid_carbons} carbons in fatty acid chain"
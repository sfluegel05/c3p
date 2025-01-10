"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside
A D-galactoside having beta-configuration at its anomeric centre
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_beta_D_galactoside, reason)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic pyranose ring pattern
    pyranose_pattern = Chem.MolFromSmarts("[OX2][C][C][C][C][C]1")
    if not mol.HasSubstructMatch(pyranose_pattern):
        return False, "No pyranose ring found"

    # Find all potential pyranose rings
    matches = mol.GetSubstructMatches(pyranose_pattern)
    
    for match in matches:
        ring_O = mol.GetAtomWithIdx(match[0])
        C1 = mol.GetAtomWithIdx(match[1])  # anomeric carbon
        C2 = mol.GetAtomWithIdx(match[2])
        C3 = mol.GetAtomWithIdx(match[3])
        C4 = mol.GetAtomWithIdx(match[4])
        C5 = mol.GetAtomWithIdx(match[5])
        
        # Check if this is potentially a galactose ring
        # 1. Each carbon should have one oxygen
        carbons = [C1, C2, C3, C4, C5]
        if not all(sum(1 for n in c.GetNeighbors() if n.GetAtomicNum() == 8) >= 1 for c in carbons):
            continue
            
        # 2. Check for required stereocenters
        if any(c.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for c in carbons):
            continue
            
        # 3. Check C4 has axial hydroxyl (characteristic of galactose)
        c4_neighbors = [n for n in C4.GetNeighbors() if n.GetAtomicNum() == 8 and n != ring_O]
        if not c4_neighbors:
            continue
            
        # 4. Check for beta configuration at anomeric center
        c1_neighbors = C1.GetNeighbors()
        non_ring_O = [n for n in c1_neighbors if n.GetAtomicNum() == 8 and n != ring_O]
        if not non_ring_O:
            continue
            
        # Check if the non-ring oxygen is connected to something (glycosidic bond)
        glycosidic_O = non_ring_O[0]
        if len(glycosidic_O.GetNeighbors()) != 2:
            continue
            
        # Verify D-galactose stereochemistry
        # C1(@H) - beta
        # C2(@H), C3(@H) - equatorial OH
        # C4(@@H) - axial OH
        # C5(@H) - maintains ring shape
        if (C1.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and
            C2.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and
            C3.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and
            C4.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW):
            
            # Check for CH2OH at C5
            c5_neighbors = [n for n in C5.GetNeighbors() if n.GetAtomicNum() == 6]
            if any(sum(1 for nn in n.GetNeighbors() if nn.GetAtomicNum() == 8) == 1 
                  for n in c5_neighbors):
                return True, "Contains beta-D-galactoside moiety with correct stereochemistry"
    
    return False, "No beta-D-galactoside configuration found"
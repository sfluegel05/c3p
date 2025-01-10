"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid 
with a hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a tetradecanoate ester group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for tetradecanoate ester
    # More specific pattern that ensures:
    # - Exactly 13 carbons in chain before ester group
    # - No branching or substitutions on the chain
    # - Proper ester linkage
    pattern = """
        [CH3X4:1]                     # Terminal methyl
        [CH2X4:2]                     # 12 methylene groups
        [CH2X4:3]
        [CH2X4:4]
        [CH2X4:5]
        [CH2X4:6]
        [CH2X4:7]
        [CH2X4:8]
        [CH2X4:9]
        [CH2X4:10]
        [CH2X4:11]
        [CH2X4:12]
        [CH2X4:13]
        [CX3:14](=[OX1:15])          # Carbonyl group
        [OX2:16]                      # Ester oxygen
        [!O&!S:17]                    # Connected to carbon (not O or S)
    """
    pattern = "".join(pattern.split())  # Remove whitespace
    
    tetradecanoate_pattern = Chem.MolFromSmarts(pattern)
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester group found"
        
    matches = mol.GetSubstructMatches(tetradecanoate_pattern)
    
    for match in matches:
        # Get the matched atoms
        chain_atoms = match[0:13]  # First 13 carbons
        carbonyl_carbon = match[13]
        carbonyl_oxygen = match[14]
        ester_oxygen = match[15]
        
        # Verify the chain is linear (no rings)
        ring_info = mol.GetRingInfo()
        if any(ring_info.IsAtomInRing(atom_idx) for atom_idx in chain_atoms):
            continue
            
        # Verify no branching on the chain
        has_branches = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count carbon neighbors
            carbon_neighbors = len([n for n in atom.GetNeighbors() 
                                 if n.GetAtomicNum() == 6 and 
                                 n.GetIdx() not in chain_atoms])
            if carbon_neighbors > 0:  # Any additional carbon means branching
                has_branches = True
                break
                
        if has_branches:
            continue
            
        # Verify each carbon in chain has correct number of hydrogens
        has_wrong_hydrogens = False
        for i, atom_idx in enumerate(chain_atoms):
            atom = mol.GetAtomWithIdx(atom_idx)
            if i == 0:  # Terminal methyl
                if atom.GetTotalNumHs() != 3:
                    has_wrong_hydrogens = True
                    break
            else:  # Methylene groups
                if atom.GetTotalNumHs() != 2:
                    has_wrong_hydrogens = True
                    break
                    
        if has_wrong_hydrogens:
            continue
            
        # Verify ester linkage
        ester_atom = mol.GetAtomWithIdx(ester_oxygen)
        if not (len(ester_atom.GetNeighbors()) == 2 and
                all(n.GetAtomicNum() == 6 for n in ester_atom.GetNeighbors())):
            continue
            
        # If we get here, we've found a valid tetradecanoate ester group
        return True, "Contains tetradecanoate (myristoyl) ester group"
        
    return False, "No valid tetradecanoate ester group found"
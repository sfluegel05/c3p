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

    # Pattern for beta-D-galactopyranose core with correct stereochemistry
    # C1 (anomeric) - beta configuration [@H]
    # C2 - equatorial OH [@H]
    # C3 - equatorial OH [@H]
    # C4 - axial OH [@@H]
    # C5 - [H] configuration maintains ring shape
    # C6 - CH2O group
    galactose_pattern = Chem.MolFromSmarts('''
        [OX2][C@H]1[C@H]([OX2])[C@H]([OX2])[C@@H]([OX2])[C@H]([C][OX2])O1
    ''')
    
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactopyranose core found"
    
    # Get matches for the galactose pattern
    matches = mol.GetSubstructMatches(galactose_pattern)
    
    for match in matches:
        # Get the anomeric carbon (C1)
        c1_idx = match[1]
        c1 = mol.GetAtomWithIdx(c1_idx)
        
        # Check for beta configuration at anomeric center
        # In beta configuration, the substituent at C1 should be equatorial
        o_neighbors = [n for n in c1.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(o_neighbors) != 2:
            continue
            
        # One oxygen should be ring oxygen, other should be glycosidic/substituent
        ring_o = None
        subst_o = None
        for o in o_neighbors:
            if len([n for n in o.GetNeighbors() if n.GetAtomicNum() == 6]) == 2:
                ring_o = o
            else:
                subst_o = o
                
        if ring_o is None or subst_o is None:
            continue
            
        # Check C4 hydroxyl is axial (characteristic of galactose)
        c4_idx = match[4]
        c4 = mol.GetAtomWithIdx(c4_idx)
        c4_neighbors = [n for n in c4.GetNeighbors() if n.GetAtomicNum() == 8]
        if not any(n.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW for n in c4_neighbors):
            continue

        # Verify all required stereocenters are present
        required_stereocenters = [c1_idx, match[2], match[3], c4_idx]
        if all(mol.GetAtomWithIdx(idx).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED 
               for idx in required_stereocenters):
            return True, "Contains beta-D-galactoside moiety with correct stereochemistry"
    
    return False, "No beta-D-galactoside configuration found"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        # Methyl beta-D-galactoside
        "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",
        # beta-D-Galp-(1->3)-D-Xylp
        "O([C@@H]1O[C@@H]([C@H](O)[C@H](O)[C@H]1O)CO)[C@@H]2[C@@H](O)C(OC[C@H]2O)O",
        # Not a galactoside (glucose)
        "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)C1=O"
    ]
    
    for smi in examples:
        result, reason = is_beta_D_galactoside(smi)
        print(f"SMILES: {smi}")
        print(f"Is beta-D-galactoside: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()
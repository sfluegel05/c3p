"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule contains an aldehyde group based on its SMILES string.
    An aldehyde has structure R-C(=O)H where:
    - One carbon is double bonded to oxygen (carbonyl)
    - That carbon has exactly one hydrogen
    - That carbon has exactly one other connection (to R group)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains aldehyde group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)
    
    # SMARTS pattern for aldehyde: [C;H1](=O)[#6,#7,#8,#9,#15,#16,#17,#35,#53]
    # Meaning:
    # - Carbon with exactly one hydrogen
    # - Double bonded to oxygen
    # - Single bonded to any reasonable atom that could be part of R group
    aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)[#6,#7,#8,#9,#15,#16,#17,#35,#53]")
    
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    # Get all matches
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # For each match, verify it's a true aldehyde
    for match in matches:
        carbon_idx = match[0]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        
        # Count number of hydrogens
        num_h = sum(1 for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        if num_h != 1:
            continue
            
        # Count total number of bonds
        num_bonds = carbon.GetTotalNumHs() + len([b for b in carbon.GetBonds()])
        if num_bonds != 3:  # Should have exactly 3 bonds (1 double to O, 1 single to H, 1 single to R)
            continue
            
        # Verify one oxygen double bond
        double_bonded_o = False
        for bond in carbon.GetBonds():
            if (bond.GetBondType() == Chem.BondType.DOUBLE and 
                bond.GetOtherAtom(carbon).GetAtomicNum() == 8):
                double_bonded_o = True
                break
                
        if double_bonded_o:
            return True, "Contains aldehyde group R-C(=O)H"
            
    return False, "No valid aldehyde group found"
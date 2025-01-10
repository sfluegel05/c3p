"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone (alpha,beta-unsaturated ketone)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule contains an enone group based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with general formula R(1)R(2)C=CR(3)-C(=O)R(4)
    where R(4) ≠ H, and the C=O is conjugated to a C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an enone group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for enone group
    # [C;!$(C=O)]: carbon not part of a carbonyl
    # [C;!$(C=O)]: another carbon not part of a carbonyl
    # =O: oxygen double bonded to carbon
    # [#6;!H0]: carbon with at least one hydrogen
    enone_pattern = Chem.MolFromSmarts("[C;!$(C=O)]=[C;!$(C=O)]-C(=O)[#6;!H0]")
    
    # Alternative pattern for cyclic enones and other variations
    enone_pattern2 = Chem.MolFromSmarts("[C;!$(C=O)]=[C;!$(C=O)]-C(=O)-[#6]")
    
    # Check for matches
    matches = mol.GetSubstructMatches(enone_pattern)
    matches2 = mol.GetSubstructMatches(enone_pattern2)
    
    all_matches = set(matches + matches2)
    
    if not all_matches:
        return False, "No enone group found"
    
    # Verify conjugation by checking bond types
    for match in all_matches:
        # Get the relevant atoms
        c1, c2, c3, _ = match
        
        # Check bond types
        bond1 = mol.GetBondBetweenAtoms(c1, c2)
        bond2 = mol.GetBondBetweenAtoms(c2, c3)
        
        # Verify double bond between C1-C2 and single bond between C2-C3
        if (bond1.GetBondType() == Chem.BondType.DOUBLE and 
            bond2.GetBondType() == Chem.BondType.SINGLE):
            
            # Get the carbonyl carbon
            carbonyl_carbon = mol.GetAtomWithIdx(c3)
            
            # Verify it's a ketone (not an aldehyde)
            if carbonyl_carbon.GetTotalNumHs() == 0:
                return True, "Contains alpha,beta-unsaturated ketone (enone) group"
    
    return False, "Structure has correct atoms but incorrect conjugation pattern"
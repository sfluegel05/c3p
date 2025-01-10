"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is obtained by condensation of myristic acid 
    (tetradecanoic acid) with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 14-carbon chain ending with a carboxyl group for tetradecanoic acid
    tetradecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)[OH1]")
    
    # Ester linkage pattern: an ester bond -COO-R'
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[*]")
    
    # Check for tetradecanoic acid moiety
    if not mol.HasSubstructMatch(tetradecanoic_acid_pattern):
        return False, "No tetradecanoic acid moiety found"
    
    # Get matches for the tetradecanoic acid moiety
    myristic_matches = mol.GetSubstructMatches(tetradecanoic_acid_pattern)

    # Check for correct ester linkage attached to a carboxyl group of tetradecanoic acid
    for match in myristic_matches:
        terminal_carbon_idx = match[-2]  # The carbon that is double-bonded to the oxygen in CO
        carboxyl_oxygen_idx = match[-1]  # The hydroxyl oxygen
        
        # Check if carboxyl carbon is involved in an ester linkage
        ester_found = False
        for atom in mol.GetAtomWithIdx(terminal_carbon_idx).GetNeighbors():
            if atom.GetIdx() == carboxyl_oxygen_idx:
                continue
            bond = mol.GetBondBetweenAtoms(terminal_carbon_idx, atom.GetIdx())
            if bond.GetBondType() == Chem.BondType.SINGLE and mol.GetAtomWithIdx(atom.GetIdx()).GetAtomicNum() == 8:
                if mol.HasSubstructMatch(ester_linkage_pattern):
                    ester_found = True
                    break
        
        if ester_found:
            return True, "Contains tetradecanoic acid moiety with ester linkage"

    return False, "Ester linkage not properly connected to tetradecanoic acid"
"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester of tetradecanoic acid (myristic acid) with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the tetradecanoyl group attached to an ester.
    # [CH2] part at the end allows for branched chains.
    tetradecanoyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[OX2]-[#6]")
    if not mol.HasSubstructMatch(tetradecanoyl_ester_pattern):
        return False, "No tetradecanoate ester linkage found"

    # Get the substructure matches of the ester bond.
    ester_matches = mol.GetSubstructMatches(tetradecanoyl_ester_pattern)
    
    # Check if the carbonyl carbon of the ester is directly attached to a tetradecanoyl group
    for match in ester_matches:
        carbonyl_carbon_index = match[0]

        # get neighboring carbons
        neighbors = mol.GetAtomWithIdx(carbonyl_carbon_index).GetNeighbors()
        
        # Find the neighbor which is a carbon with 3 bonds (not the ester oxygen)
        acyl_carbon_neighbor = None
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalDegree() == 3:
                acyl_carbon_neighbor = neighbor
                break
        
        if acyl_carbon_neighbor is None:
             return False, "No tetradecanoyl group attached to the carbon of an ester group found"
        
        
        # Define a SMARTS pattern for a tetradecanoyl group chain starting at the carbonyl carbon
        # The pattern requires 13 more carbons
        tetradecanoyl_chain_pattern = Chem.MolFromSmarts(f"[CX3]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") 

        
        # Check if tetradecanoyl carbon has a long chain
        if not mol.HasSubstructMatch(tetradecanoyl_chain_pattern,  useIndices=True, atomMap={0:acyl_carbon_neighbor.GetIdx()}):
            return False, "Not a tetradecanoyl chain attached to the carbon of an ester group"
        
        # Count the number of carbons directly attached to the carbonyl of the ester.
        carbon_count=0
        
        for atom in mol.GetSubstructMatch(tetradecanoyl_chain_pattern, useIndices=True, atomMap={0:acyl_carbon_neighbor.GetIdx()}):
            if atom.GetAtomicNum()==6:
                carbon_count += 1
        
        if carbon_count != 13:
              return False, "Not a 13 carbon chain on the tetradecanoyl"

    return True, "Contains tetradecanoyl group in an ester linkage"
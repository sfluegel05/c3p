"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester consists of an ester linkage formed between a fatty acid and a fatty alcohol.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Find ester linkage(s)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Verify carbon chain lengths (fatty acid/alcohol characteristics)
    # Specific length requirements are not well-defined universally, but they are typically long for wax esters
    for match in ester_matches:
        # Check carbon chain before and after the ester group
        atom_idx_o = match[2]  # Index of the oxygen atom involved in the ester
        atom_idx_c = match[0]  # Index of the carbon atom double-bonded to oxygen
        chain_before = Chem.GetShortestPath(mol, atom_idx_c, 0)  # Assumes chain starts from this carbon
        chain_after = Chem.GetShortestPath(mol, atom_idx_o, mol.GetNumAtoms() - 1)  # Assumes chain ends at the last atom

        # Estimate the length, ensuring more than a minimal length typical for fatty acids/alcohols
        if len(chain_before) < 8 or len(chain_after) < 8:
            return False, f"Chains are too short to be characteristic of a wax ester (chains found: {len(chain_before)} and {len(chain_after)})"

    return True, "Contains ester linkage with suitable long chains typical for wax esters"
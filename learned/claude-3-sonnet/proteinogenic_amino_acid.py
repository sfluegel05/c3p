"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:33709 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is an alpha-amino acid
    patt = Chem.MolFromSmarts('N[C@@H](C(=O)O)C')
    if not mol.HasSubstructMatch(patt):
        return False, "Not an alpha-amino acid"
    
    # Check the side chain
    patt = Chem.MolFromSmarts('N[C@@H](C(=O)O)[C@@H]')
    if not mol.HasSubstructMatch(patt):
        return False, "No side chain attached"
    
    # Get the side chain
    sidechains = AllChem.GetMolFragsWithSPorts(mol, AllChem.SmilesPatts.smarts_pattern['sidechain'])
    if not sidechains:
        return False, "Could not determine side chain"
    
    # Check if the side chain matches a proteinogenic amino acid
    for chain in sidechains:
        smarts = AllChem.GetSmilesPatternFromMol(chain)
        if smarts in AllChem.SmartsPatts.proteinogenic_amino_acids.values():
            return True, "Matches a proteinogenic amino acid side chain"
    
    return False, "Side chain does not match a proteinogenic amino acid"
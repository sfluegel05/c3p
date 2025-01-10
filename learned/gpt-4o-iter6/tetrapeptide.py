"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide must contain four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for the generic amide bond pattern and count them
    amide_pattern = Chem.MolFromSmarts("N-[C;$(C(=O))]-[!#1]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 3:
        return False, f"Contains {len(amide_matches)} peptide bonds, expected 3"
    
    # Determine the number of amino acid residues (look for repeating CN backbone units)
    residue_count = 0
    for match in amide_matches:
        carbon_idx, nitrogen_idx = match[1], match[0]
        carbon = mol.GetAtomWithIdx(carbon_idx)
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Verify if these match the typical amino acid backbone (considering variable side chains)
        if carbon.GetDegree() >= 3 and nitrogen.GetDegree() >= 2:
            residue_count += 1
    
    if residue_count != 4:
        return False, f"Peptide linkages found, but amino acid count is {residue_count}, expected 4"

    return True, "Contains four amino-acid residues connected by peptide linkages"
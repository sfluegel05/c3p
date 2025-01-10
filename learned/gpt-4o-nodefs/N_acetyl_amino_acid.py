"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is an N-acetyl-amino acid, False otherwise
                             along with a reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetyl group ('CC(=O)N') attached directly to nitrogen
    acetyl_pattern = Chem.MolFromSmarts('[#6]-C(=O)-N')
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    if not acetyl_matches:
        return False, "No acetyl group attached to nitrogen found"
    
    # Look for carboxylic acid group attached in the amino acid structure
    carboxyl_pattern = Chem.MolFromSmarts(Chem.MolToSmarts(Chem.MolFromSmiles('C(=O)O')))
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    # Ensure the characteristic backbone of amino acid is present: N-CH-C(=O)O structure
    # specifically looking for N-[C@@H]-C(=O)O
    amino_acid_pattern = Chem.MolFromSmarts('[NX3][C@@H]-C(=O)O')
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid backbone identified"
    
    # Final verification: ensure acetyl and amino acid structures are connected
    for acetyl_match in acetyl_matches:
        acetyl_nitrogen_idx = next((idx for idx in acetyl_match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7), None)
        if acetyl_nitrogen_idx is not None:
            for nei in mol.GetAtomWithIdx(acetyl_nitrogen_idx).GetNeighbors():
                if nei.GetIdx() in {idx for subset in amino_acid_matches for idx in subset}:
                    return True, "Contains acetyl group correctly bound to nitrogen of an amino acid"
                    
    return False, "Acetyl group not correctly bound to the nitrogen of an amino acid"
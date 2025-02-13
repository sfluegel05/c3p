"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify inositol ring (6-membered carbocycl with stereochemistry and multiple hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("C1([C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)O)C(O)C(O)C(O)C(O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"
    
    # Identify phosphodiester bridge
    phospho_pattern = Chem.MolFromSmarts("O=P(O)(OC)OC")
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "No phosphodiester bridge found"
    
    # Identify ceramide moiety (amide bond plus long chain)
    ceramide_pattern = Chem.MolFromSmarts("N[C@@H](CO)C(=O)C")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # Additional check for long alkyl chain in ceramide
    alkane_chain_length = sum([1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and any(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors())])
    if alkane_chain_length < 12:
        return False, f"Alkyl chain length in ceramide too short: {alkane_chain_length}"
    
    return True, "Pattern matches inositol phosphoceramide structure"
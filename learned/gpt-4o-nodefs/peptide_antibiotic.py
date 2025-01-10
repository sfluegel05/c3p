"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide Antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic typically contains multiple amide bonds, cyclic structures,
    and functional groups indicative of bioactive peptides.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a potential peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count amide bonds (R-C(=O)-NR) and require more than 5
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 5:
        return False, "Insufficient amide (peptide) bonds"

    # Check for complex cyclic structures - ensure there's more than one ring
    num_rings = rdmolops.GetSSSR(mol)
    if num_rings < 1:
        return False, "No cyclic structures detected"

    # Molecular features like specific ring systems or unusual amino acids could be added here
    # For example, adding detection for Thiazole ring systems which are common in some peptide antibiotics
    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")
    if not mol.HasSubstructMatch(thiazole_pattern):
        return False, "No thiazole rings detected; common in some peptide antibiotics"

    # Molecular weight check - peptide antibiotics are generally large 
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for complex peptide antibiotic"

    return True, "Structural features typical of a peptide antibiotic found"
"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Common nucleobase scaffolds
    pyrimidine_patterns = [
        Chem.MolFromSmarts("c1ncnc1"),  # Basic pyrimidine structure
        Chem.MolFromSmarts("c1cncnc1"),  # Extended variants of pyrimidine
        Chem.MolFromSmarts("c1ncnn1")  # Other nitrogen substitutions
    ]
    purine_patterns = [
        Chem.MolFromSmarts("c1ncnc2[nH]ncnc12"),  # Basic purine structure
        Chem.MolFromSmarts("c1ncnc2ncnc[nH]12"),  # Extended variants of purine
        Chem.MolFromSmarts("c1ncnc2ncnn12"),  # Other nitrogen substitutions
    ]
    
    # Check for pyrimidine or purine-like structures
    for pyrimidine in pyrimidine_patterns:
        if mol.HasSubstructMatch(pyrimidine):
            return True, "Contains pyrimidine-like heterocycle"
            
    for purine in purine_patterns:
        if mol.HasSubstructMatch(purine):
            return True, "Contains purine-like heterocycle"
    
    # Look for key functional groups common to nucleobases
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # Carbonyl group
    amino_pattern = Chem.MolFromSmarts("[NX3][H]")  # Amino group
    
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Check for required functional group combinations
    if len(carbonyl_matches) > 0 and len(amino_matches) > 0:
        return True, "Contains both carbonyl and amino groups typical of nucleobases"
    
    # Analyze nitrogen content and heteroatom features
    n_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    hetero_atm_count = rdMolDescriptors.CalcNumHeteroatoms(mol)

    if n_atoms >= 2 and hetero_atm_count >= 3:
        return True, "Complex nitrogen and heteroatom structure suggesting a modified nucleobase"
    
    return False, "Does not match features of known nucleobase analogues"
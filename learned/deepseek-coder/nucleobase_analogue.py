"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
Definition: A molecule that can substitute for a normal nucleobase in nucleic acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue must resemble a natural nucleobase (adenine, guanine, cytosine, thymine, uracil)
    but with modifications such as substitutions, additions, or alterations.

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

    # Define more specific SMARTS patterns for nucleobase cores
    purine_core = Chem.MolFromSmarts("[nH]1[c,n][c,n][c,n][c,n][c,n]1")  # Purine core
    pyrimidine_core = Chem.MolFromSmarts("[nH]1[c,n][c,n][c,n]1")  # Pyrimidine core

    # Check if the molecule has a nucleobase-like core structure
    has_core = (mol.HasSubstructMatch(purine_core) or 
                mol.HasSubstructMatch(pyrimidine_core))
    
    if not has_core:
        return False, "No nucleobase-like core structure found"

    # Check for common nucleobase functional groups
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H2,H1]")  # Amide group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # Carbonyl group
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")  # Amine group

    has_functional_groups = (mol.HasSubstructMatch(amide_pattern) or
                            mol.HasSubstructMatch(carbonyl_pattern) or
                            mol.HasSubstructMatch(amine_pattern))
    
    # Check for common modifications in nucleobase analogues
    modification_patterns = [
        Chem.MolFromSmarts("[#6]=[#6]"),  # Double bonds (etheno bridges)
        Chem.MolFromSmarts("[#7]-[#7]"),  # Azide or diazo groups
        Chem.MolFromSmarts("[#6]-[#7]=[#7]"),  # Triazenes
        Chem.MolFromSmarts("[#6]-[#8]-[#6]"),  # Ether linkages
        Chem.MolFromSmarts("[#6]-[#16]-[#6]"),  # Thioether linkages
        Chem.MolFromSmarts("[#6]-[#17]"),  # Halogen substitutions
        Chem.MolFromSmarts("[#6]-[#9]"),  # Fluorine substitutions
        Chem.MolFromSmarts("[#6]-[#35]"),  # Bromine substitutions
        Chem.MolFromSmarts("[#6]-[#53]")  # Iodine substitutions
    ]
    
    has_modifications = any(mol.HasSubstructMatch(patt) for patt in modification_patterns)
    
    # Check molecular weight (nucleobase analogues are typically small to medium molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:  # Increased threshold
        return False, "Molecular weight too high for nucleobase analogue"

    # If it has the core structure and either functional groups or modifications, it's likely an analogue
    if has_core and (has_functional_groups or has_modifications):
        return True, "Resembles a nucleobase analogue with modifications"
    
    return False, "Does not resemble a nucleobase analogue"
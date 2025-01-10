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

    # Define more flexible SMARTS patterns for nucleobase cores
    purine_like = Chem.MolFromSmarts("[n]1[c,n][c,n][c,n][c,n][c,n]1")  # Purine-like 6+5 ring system
    pyrimidine_like = Chem.MolFromSmarts("[n]1[c,n][c,n][c,n]1")  # Pyrimidine-like 6-ring system
    imidazole_like = Chem.MolFromSmarts("[n]1[c,n][c,n]1")  # Imidazole-like 5-ring system

    # Check if the molecule has a nucleobase-like core structure
    has_core = (mol.HasSubstructMatch(purine_like) or 
                mol.HasSubstructMatch(pyrimidine_like) or
                mol.HasSubstructMatch(imidazole_like))
    
    if not has_core:
        return False, "No nucleobase-like core structure found"

    # Check for common nucleobase features
    # Nitrogen atoms in the core (at least 2)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 2:
        return False, "Insufficient nitrogen atoms for nucleobase analogue"

    # Check for common functional groups in nucleobases
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3]")
    
    has_functional_groups = (mol.HasSubstructMatch(amine_pattern) or
                            mol.HasSubstructMatch(carbonyl_pattern) or
                            mol.HasSubstructMatch(imine_pattern))
    
    if not has_functional_groups:
        return False, "No characteristic nucleobase functional groups found"

    # Check molecular weight (nucleobase analogues are typically small molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for nucleobase analogue"

    # Check for common modifications in nucleobase analogues
    modification_patterns = [
        Chem.MolFromSmarts("[#6]=[#6]"),  # Double bonds (etheno bridges)
        Chem.MolFromSmarts("[#7]-[#7]"),  # Azide or diazo groups
        Chem.MolFromSmarts("[#6]-[#7]=[#7]"),  # Triazenes
        Chem.MolFromSmarts("[#6]-[#8]-[#6]"),  # Ether linkages
        Chem.MolFromSmarts("[#6]-[#16]-[#6]")  # Thioether linkages
    ]
    
    has_modifications = any(mol.HasSubstructMatch(patt) for patt in modification_patterns)
    
    # If it has the core structure and either functional groups or modifications, it's likely an analogue
    if has_core and (has_functional_groups or has_modifications):
        return True, "Resembles a nucleobase analogue with modifications"
    
    return False, "Does not resemble a nucleobase analogue"
"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: CHEBI:33587 organochlorine compound
An organochlorine compound is a compound containing at least one carbon-chlorine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound is a compound containing at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carbon-chlorine bonds
    pattern = Chem.MolFromSmarts("[Cl;X1]-[#6]")
    has_carbon_chlorine_bond = mol.HasSubstructMatch(pattern)
    
    if not has_carbon_chlorine_bond:
        return False, "No carbon-chlorine bonds found"
    
    # Check for common non-organochlorine substructures (to avoid false positives)
    excluded_patterns = [
        Chem.MolFromSmarts("[Cl-]"),  # Inorganic chloride anion
        Chem.MolFromSmarts("[Cl+]"),  # Chlorine cation
        Chem.MolFromSmarts("[Cl]~[Cl]"),  # Chlorine gas
        Chem.MolFromSmarts("[Cl]~[#8]~[#8]"),  # Hypochlorite
        Chem.MolFromSmarts("[Cl]~[#6](=[#8])~[#8]"),  # Chloroformate
        Chem.MolFromSmarts("[Cl]~[#6](=[#8])[#8]"),  # Chloroacetate
    ]
    
    for pattern in excluded_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Excluded non-organochlorine substructure found"
    
    # Additional checks for common organochlorine functional groups
    organochlorine_patterns = [
        Chem.MolFromSmarts("[Cl]-[#6]-[#6]"),  # Alkyl chloride
        Chem.MolFromSmarts("[Cl]-[#6]-[#6]-[#8]"),  # Chloroalkyl ether
        Chem.MolFromSmarts("[Cl]-[#6]-[#6]=[#8]"),  # Chloroalkyl ketone
        Chem.MolFromSmarts("[Cl]-[#6]=[#6]"),  # Vinyl chloride
        Chem.MolFromSmarts("[Cl]-[#6]~[#6]~[#6]~[#6]~[#6]"),  # Chlorinated cyclic alkane
        Chem.MolFromSmarts("[Cl]-[#6]~[#6]~[#6]~[#6]~[#6]~[#6]"),  # Chlorinated aromatic ring
    ]
    
    has_organochlorine_group = False
    for pattern in organochlorine_patterns:
        if mol.HasSubstructMatch(pattern):
            has_organochlorine_group = True
            break
    
    if not has_organochlorine_group:
        return False, "No common organochlorine functional groups found"
    
    # Additional check for molecular weight (organochlorines typically > 100 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for an organochlorine compound"
    
    return True, "Molecule contains at least one carbon-chlorine bond and common organochlorine functional groups"
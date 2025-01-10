"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: catecholamine
Definition: 4-(2-Aminoethyl)pyrocatechol and derivatives formed by substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to neutral form if possible (handle salt forms)
    mol = Chem.RemoveHs(mol)  # Remove explicit hydrogens
    
    # Look for catechol (benzene with two adjacent OH groups) pattern
    # Allow for additional OH groups and other substituents
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cccc1")  # Basic catechol
    catechol_pattern2 = Chem.MolFromSmarts("c1c(O)cc(O)ccc1")  # Alternative OH positions
    catechol_pattern3 = Chem.MolFromSmarts("c1c(O)c(O)cc(O)c1")  # Trihydroxy pattern
    
    has_catechol = any([mol.HasSubstructMatch(pat) for pat in 
                       [catechol_pattern, catechol_pattern2, catechol_pattern3]])
    
    if not has_catechol:
        return False, "No catechol (adjacent hydroxyl groups on benzene) found"
    
    # Look for ethylamine chain (-CH2-CH2-N) attached to the ring
    # Allow for substituted amines and various attachment points
    ethylamine_patterns = [
        Chem.MolFromSmarts("c-CCN"),  # Basic ethylamine
        Chem.MolFromSmarts("c-CC[NH2]"),  # Primary amine
        Chem.MolFromSmarts("c-CC[NH][CH3]"),  # Secondary amine (methylated)
        Chem.MolFromSmarts("c-CC[NH]C"),  # Secondary amine (general)
        Chem.MolFromSmarts("c-CCN(C)C"),  # Tertiary amine
        Chem.MolFromSmarts("c-CC[NH+]"),  # Protonated amine
        Chem.MolFromSmarts("c-C[CH]N"),  # Allow for substitution on carbon
        Chem.MolFromSmarts("c-C[CH](O)N"),  # Beta-hydroxyl variant
    ]
    
    has_ethylamine = any([mol.HasSubstructMatch(pat) for pat in ethylamine_patterns])
    
    if not has_ethylamine:
        return False, "No ethylamine chain found"
    
    # Additional check to ensure the catechol and ethylamine are properly connected
    # by looking for the complete pattern
    complete_pattern = Chem.MolFromSmarts("c1c(O)c(O)ccc1CCN")  # Basic pattern
    complete_pattern2 = Chem.MolFromSmarts("c1c(O)c(O)cc(CCN)c1")  # Alternative attachment
    
    has_complete = any([mol.HasSubstructMatch(pat) for pat in 
                       [complete_pattern, complete_pattern2]])
    
    if not has_complete:
        return False, "Catechol and ethylamine groups not properly connected"
        
    # Check for reasonable molecular size (to exclude very large molecules)
    if mol.GetNumAtoms() > 50:
        return False, "Molecule too large to be a typical catecholamine"
    
    return True, "Contains catechol moiety with ethylamine chain"
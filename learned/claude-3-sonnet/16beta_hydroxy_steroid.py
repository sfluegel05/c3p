"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:16beta-hydroxy steroid
A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for stereochemistry
    mol = Chem.AddHs(mol)
    
    try:
        # More flexible steroid core pattern
        # Four connected rings with some flexibility in saturation
        steroid_core = Chem.MolFromSmarts("C1~C~C~C2~C~C~C3~C~C~C4~C~C~C~C4~C3~C2~1")
        if steroid_core is None:
            return False, "Invalid steroid core SMARTS pattern"
            
        if not mol.HasSubstructMatch(steroid_core):
            return False, "No steroid core structure found"

        # Pattern for 16-beta-hydroxy group
        # More specific pattern focusing on position 16 and beta configuration
        # [C@@H] indicates specific stereochemistry
        beta_oh_16 = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]2~[C]~[C]~[C]3~[C]~[C]~[C]4~[C]~[C]~[C@@H](O)~[C]4~[C]3~[C]2~1")
        if beta_oh_16 is None:
            return False, "Invalid 16-beta-hydroxy SMARTS pattern"
            
        if not mol.HasSubstructMatch(beta_oh_16, useChirality=True):
            return False, "No 16-beta-hydroxy group found or wrong stereochemistry"

        # Basic molecular property checks
        mol_weight = Chem.Descriptors.ExactMolWt(mol)
        if mol_weight < 250 or mol_weight > 1000:
            return False, f"Molecular weight {mol_weight} outside typical steroid range (250-1000)"

        # Count carbons (steroids typically have 17+ carbons)
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count < 17:
            return False, f"Too few carbons ({carbon_count}) for a steroid structure"

        # Count oxygens (need at least one for the hydroxy group)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        if oxygen_count < 1:
            return False, "No oxygen atoms found"

        # Look for specific OH group
        hydroxy_pattern = Chem.MolFromSmarts("[OH]")
        if hydroxy_pattern is None:
            return False, "Invalid hydroxy SMARTS pattern"
            
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "No hydroxy group found"

        return True, "Contains steroid core with 16-beta-hydroxy group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"
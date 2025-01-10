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
        # Define multiple steroid core patterns to catch different variants
        steroid_patterns = [
            # Basic steroid core (more flexible)
            "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1",
            # Alternative pattern with more specific bond types
            "C~1~C~C~C~2~C~C~C~3~C~C~C~4~C~C~C(O)~C~4~C~3~C~2~1",
            # Pattern for estrane derivatives
            "c1cc2C~C~C~3~C~C~C~4~C~C~C(O)~C~4~C~3~C~c2cc1"
        ]
        
        has_steroid_core = False
        for pattern in steroid_patterns:
            steroid_core = Chem.MolFromSmarts(pattern)
            if steroid_core is not None and mol.HasSubstructMatch(steroid_core):
                has_steroid_core = True
                break
                
        if not has_steroid_core:
            return False, "No steroid core structure found"

        # Pattern specifically for 16-beta-hydroxy
        # This pattern looks for the D ring with a beta-OH at position 16
        # The [C@@H] ensures beta stereochemistry
        beta_oh_patterns = [
            # Pattern for saturated D ring with 16β-OH
            "[C]~1~[C]~[C]~[C]~2~[C]~[C]~[C]~3~[C@@H](O)~[C]~[C]~[C]~3~[C]~2~1",
            # Alternative pattern with more explicit stereochemistry
            "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[C@@H](O)~[#6]~[#6]~[#6]~3~[#6]~2~1"
        ]
        
        has_16beta_oh = False
        for pattern in beta_oh_patterns:
            beta_oh = Chem.MolFromSmarts(pattern)
            if beta_oh is not None and mol.HasSubstructMatch(beta_oh, useChirality=True):
                has_16beta_oh = True
                break
                
        if not has_16beta_oh:
            return False, "No 16-beta-hydroxy group found or wrong stereochemistry"

        # Basic molecular property checks
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count < 17:
            return False, f"Too few carbons ({carbon_count}) for a steroid structure"

        # Count oxygens (need at least one for the hydroxy group)
        oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        if oxygen_count < 1:
            return False, "No oxygen atoms found"

        # Additional check for hydroxyl group
        hydroxy_pattern = Chem.MolFromSmarts("[OH]")
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "No hydroxy group found"

        return True, "Contains steroid core with 16-beta-hydroxy group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"
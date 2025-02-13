"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol compounds
Definition: A chromanol with a chroman-6-ol skeleton substituted at position 2 
by a saturated or triply-unsaturated hydrocarbon chain of three isoprenoid units
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_tocol, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for chromanol core - including both saturated and unsaturated variants
    # Also matches the required hydroxyl position
    chromanol_patterns = [
        # Saturated chromanol core with hydroxyl
        "O1CCCc2c(c(O)c(C)cc2)C1",
        # Alternative pattern for variations
        "O1CCCc2c(c(O)ccc2)C1",
        # Pattern for unsaturated pyran ring
        "O1C=CCc2c(c(O)c(C)cc2)C1"
    ]
    
    core_found = False
    for pattern in chromanol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            core_found = True
            break
    
    if not core_found:
        return False, "No valid chromanol core found"
    
    # Check for proper substitution at position 2
    # This ensures the side chain is attached at the correct position
    c2_substitution = Chem.MolFromSmarts("O1[C]([CH2,CH3])[CH2,CH]Cc2c(c(O)ccc2)C1")
    if not mol.HasSubstructMatch(c2_substitution):
        return False, "Incorrect substitution at position 2"
    
    # Check for proper isoprenoid chain - both saturated and unsaturated versions
    isoprenoid_patterns = [
        # Saturated chain pattern
        "[CH3]C(C)CCC[CH](C)CCC[CH](C)CCC",
        # Unsaturated chain pattern
        "[CH3]C(C)CCC([CH3])=CCC([CH3])=CCC",
        # Alternative unsaturated pattern
        "[CH3]C(C)=CCC([CH3])=CCC([CH3])=CCC"
    ]
    
    chain_found = False
    for pattern in isoprenoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            chain_found = True
            break
    
    if not chain_found:
        return False, "Missing required isoprenoid chain structure"
    
    # Check for proper carbon count and overall structure
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 20:  # Minimum carbons for a tocol
        return False, "Insufficient carbons for tocol structure"
    
    # Check for allowed modifications without affecting classification
    modifications = []
    mod_patterns = {
        "acetate": "CC(=O)O",
        "succinate": "O=C(O)CCC(=O)O",
        "carboxyl": "C(=O)O"
    }
    
    for mod_name, pattern in mod_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            modifications.append(mod_name)
    
    # Additional check to prevent false positives with modified chains
    # Exclude molecules with hydroxyl groups in the side chain
    side_chain_oh = Chem.MolFromSmarts("O1CCCc2ccccc12CC[CH2][OH]")
    if mol.HasSubstructMatch(side_chain_oh):
        return False, "Invalid modification in side chain"
    
    base_message = "Contains chroman-6-ol core with appropriate isoprenoid chain"
    if modifications:
        return True, f"{base_message} (modified with: {', '.join(modifications)})"
    
    return True, base_message
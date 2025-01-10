"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count basic atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 10:
        return False, "Too few carbons for cannabinoid structure"
    if num_oxygens < 1:
        return False, "Cannabinoids must contain oxygen"

    # Classical cannabinoid patterns (like THC)
    classical_pattern = Chem.MolFromSmarts("[OX2H1]c1c(C)cc(CCCCC)cc1OC") # Simplified THC-like core
    
    # Endocannabinoid patterns
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO") # Ethanolamine group
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO") # Glycerol group
    long_chain = Chem.MolFromSmarts("CCCCCCCC") # At least 8 carbons in chain
    
    # Synthetic cannabinoid patterns
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1") # Indole core
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    
    # Check for characteristic functional groups
    phenol_pattern = Chem.MolFromSmarts("[OX2H1]c1ccccc1")
    ether_pattern = Chem.MolFromSmarts("[OX2](C)C")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    
    # Look for key structural features
    has_classical = mol.HasSubstructMatch(classical_pattern)
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_long_chain = mol.HasSubstructMatch(long_chain)
    has_indole = mol.HasSubstructMatch(indole_pattern)
    has_phenol = mol.HasSubstructMatch(phenol_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_benzene = mol.HasSubstructMatch(benzene_pattern)
    
    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Classification logic
    if has_classical:
        return True, "Contains classical cannabinoid core structure"
        
    if (has_ethanolamine or has_glycerol) and has_long_chain:
        if has_carbonyl:
            return True, "Matches endocannabinoid pattern with long chain and polar head group"
            
    if has_indole and has_carbonyl and ring_count >= 2:
        return True, "Matches synthetic cannabinoid pattern with indole core"
        
    # Additional checks for other cannabinoid-like structures
    if (has_phenol or has_ether) and has_long_chain:
        if ring_count >= 1 and rotatable_bonds >= 4:
            if mol_wt > 200 and mol_wt < 600:
                return True, "Contains cannabinoid-like features (rings, chains, and oxygen-containing groups)"
                
    # Check for specific structural combinations
    structural_features = sum([has_phenol, has_ether, has_carbonyl, has_benzene])
    if structural_features >= 2 and has_long_chain and ring_count >= 1:
        if 200 < mol_wt < 600:
            return True, "Contains multiple cannabinoid structural features"
            
    return False, "Does not match cannabinoid structural patterns"
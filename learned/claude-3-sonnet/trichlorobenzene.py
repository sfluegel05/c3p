"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene
    (benzene ring with exactly 3 chlorine substituents and minimal other substitution).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Quick check for molecular composition
    if rdMolDescriptors.CalcMolFormula(mol) == "C6H3Cl3":
        return True, "Molecular formula matches trichlorobenzene (C6H3Cl3)"
    
    # Find all benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"
    
    # Patterns to exclude
    dibenzofuran_pattern = Chem.MolFromSmarts("c1ccc2oc3ccccc3c2c1")
    dibenzodioxin_pattern = Chem.MolFromSmarts("c1ccc2Oc3ccccc3Oc2c1")
    
    if mol.HasSubstructMatch(dibenzofuran_pattern) or mol.HasSubstructMatch(dibenzodioxin_pattern):
        return False, "Molecule contains dibenzofuran or dibenzodioxin structure"
    
    # For each benzene ring, check its substitution pattern
    for benzene_atoms in benzene_matches:
        chlorine_count = 0
        total_substituents = 0
        has_complex_substituents = False
        
        # Check each atom in the benzene ring
        for ring_atom_idx in benzene_atoms:
            ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
            
            # Look at neighbors of ring atom
            for neighbor in ring_atom.GetNeighbors():
                if neighbor.GetIdx() not in benzene_atoms:
                    total_substituents += 1
                    
                    if neighbor.GetAtomicNum() == 17:  # Chlorine
                        chlorine_count += 1
                    elif neighbor.GetAtomicNum() not in [1, 17]:  # Not H or Cl
                        has_complex_substituents = True
                        
                    # Check if neighbor is part of a complex group
                    if len(list(neighbor.GetNeighbors())) > 1:
                        has_complex_substituents = True
        
        # Conditions for a valid trichlorobenzene:
        # 1. Exactly 3 chlorines
        # 2. Total substituents should be 3 or 4 (allowing for one additional small group)
        # 3. No complex substituents
        if (chlorine_count == 3 and 
            total_substituents <= 4 and 
            not has_complex_substituents):
            
            return True, "Found benzene ring with exactly 3 chlorine substituents and minimal other substitution"
            
    return False, "No suitable trichlorobenzene structure found"
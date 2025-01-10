"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam compounds
A lactam in which the amide bond is contained within a four-membered ring, 
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Handle salts and charged species by splitting on '.'
    fragments = smiles.split('.')
    for fragment in fragments:
        mol = Chem.MolFromSmiles(fragment)
        if mol is None:
            continue
            
        # Basic beta-lactam pattern: 4-membered ring containing N-C(=O)
        # [#7] = any nitrogen
        # ~2 = exactly 2 bonds
        # [#6] = any carbon
        # [#8] = any oxygen
        beta_lactam_pattern = Chem.MolFromSmarts('[#7]1~2[#6][#6](=[#8])[#6]1')
        
        if mol.HasSubstructMatch(beta_lactam_pattern):
            # Additional validation: confirm it's part of a 4-membered ring
            ring_info = mol.GetRingInfo()
            if any(len(ring) == 4 for ring in ring_info.AtomRings()):
                # Check for specific beta-lactam types
                patterns = [
                    ('[#7]1[#6]2[#6](=[#8])[#6]1[#16][#6]', 'penicillin/cephalosporin-type'),
                    ('[#7]1[#6]2[#6](=[#8])[#6]1[#6][#6]2', 'carbapenem-type'),
                    ('[#7]1[#6][#6](=[#8])[#6]1[#16](=[#8])(=[#8])[#8]', 'monobactam-type')
                ]
                
                for pattern, desc in patterns:
                    pat = Chem.MolFromSmarts(pattern)
                    if pat and mol.HasSubstructMatch(pat):
                        return True, f"Contains beta-lactam structure ({desc})"
                        
                return True, "Contains beta-lactam ring with amide bond"
                
    return False, "No beta-lactam ring found"
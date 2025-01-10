"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide compounds
An anilide is an aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide has a phenyl ring connected to an amide nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for N-phenyl amide pattern
    # [c] - aromatic carbon (part of phenyl ring)
    # [NX3] - nitrogen with 3 connections
    # [CX3](=[OX1]) - amide carbon with double bonded oxygen
    anilide_pattern = Chem.MolFromSmarts('[c][NX3][CX3](=[OX1])')
    
    matches = mol.GetSubstructMatches(anilide_pattern)
    if not matches:
        return False, "No N-phenyl amide group found"
        
    # For each match, verify the aromatic ring is a proper phenyl or substituted phenyl
    for match in matches:
        arom_c = match[0]  # First atom in match is the aromatic carbon
        ring = Chem.GetSymmSSSR(mol)
        
        # Check each ring that contains our aromatic carbon
        for r in ring:
            if arom_c in r:
                # Get all atoms in the ring
                ring_atoms = set(r)
                ring_size = len(ring_atoms)
                
                # Count aromatic atoms
                aromatic_count = sum(1 for atom_idx in ring_atoms 
                                   if mol.GetAtomWithIdx(atom_idx).GetIsAromatic())
                
                # Verify it's a 6-membered aromatic ring
                if ring_size == 6 and aromatic_count == 6:
                    # Verify ring atoms are C or substituted C
                    is_phenyl = True
                    for atom_idx in ring_atoms:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetAtomicNum() != 6:  # not carbon
                            is_phenyl = False
                            break
                    
                    if is_phenyl:
                        return True, "Contains N-phenyl amide group"
                        
    return False, "No proper phenyl ring attached to amide nitrogen"
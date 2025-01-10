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
    # Match both secondary and tertiary amides where N is connected to phenyl
    # [c] - aromatic carbon
    # [NX3] - nitrogen with 3 connections (includes both NH and N-substituted)
    # [CX3](=[OX1]) - amide carbon with double bonded oxygen
    # Exclude cases where amide C is bonded to N or O (to avoid urea/carbamate)
    anilide_pattern = Chem.MolFromSmarts('[c][NX3]([#1,#6])[CX3](=[OX1])[!N;!O]')
    
    matches = mol.GetSubstructMatches(anilide_pattern)
    if not matches:
        return False, "No N-phenyl amide group found"
        
    # For each match, verify the aromatic ring is a proper phenyl
    for match in matches:
        arom_c = match[0]  # First atom in match is the aromatic carbon
        ring = Chem.GetSymmSSSR(mol)
        
        # Check each ring that contains our aromatic carbon
        for r in ring:
            if arom_c in r:
                # Get all atoms in the ring
                ring_atoms = set(r)
                # Count carbons and check ring size
                ring_size = len(ring_atoms)
                carbon_count = sum(1 for atom_idx in ring_atoms 
                                 if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
                
                # Verify it's a 6-membered aromatic ring with all carbons
                if ring_size == 6 and carbon_count == 6:
                    # Check if all atoms in ring are aromatic
                    all_aromatic = all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() 
                                     for atom_idx in ring_atoms)
                    if all_aromatic:
                        return True, "Contains N-phenyl amide group"
                        
    return False, "No proper phenyl ring attached to amide nitrogen"
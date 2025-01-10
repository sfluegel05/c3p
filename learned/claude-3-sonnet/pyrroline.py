"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:51262 pyrroline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule contains a pyrroline ring based on its SMILES string.
    A pyrroline is a dihydropyrrole - a 5-membered heterocyclic ring containing 
    one nitrogen atom and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a pyrroline ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Make a copy for Kekulization
    mol_kekulized = Chem.Mol(mol)
    try:
        Chem.Kekulize(mol_kekulized, clearAromaticFlags=True)
    except:
        pass  # If kekulization fails, we'll still try with the original molecule
    
    # Multiple SMARTS patterns to catch different forms of pyrroline
    patterns = [
        # Basic pyrroline patterns with double bond in different positions
        "[NX3R5]1[#6][#6]=[#6][#6]1",  # Double bond at 3-4 position
        "[NX3R5]1[#6]=[#6][#6][#6]1",  # Double bond at 2-3 position
        "[#6]1[#6][NX3R5][#6]=[#6]1",  # Double bond at 4-5 position
        
        # Patterns for keto forms
        "[NX3R5]1[#6][#6](=[O,S])[#6][#6]1",
        "[NX3R5]1[#6](=[O,S])[#6][#6][#6]1",
        
        # Pattern for charged forms
        "[N+X4R5]1[#6][#6][#6][#6]1",
        
        # Pattern for imine forms
        "[NX2R5]=[#6R5][#6R5][#6R5][#6R5]1"
    ]

    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt) or (mol_kekulized and mol_kekulized.HasSubstructMatch(patt)):
            # Additional validation to ensure it's not fully aromatic
            rings = mol.GetRingInfo().AtomRings()
            for ring in rings:
                if len(ring) == 5:  # Check 5-membered rings
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                    if n_count == 1:  # One nitrogen
                        # Check if ring is not fully aromatic
                        aromatic_atoms = sum(1 for atom in ring_atoms if atom.GetIsAromatic())
                        if aromatic_atoms != 5:  # Not fully aromatic
                            return True, "Contains pyrroline ring (5-membered ring with N and appropriate unsaturation)"
            
            # If we found pattern but ring was aromatic
            continue
    
    return False, "No pyrroline ring found"
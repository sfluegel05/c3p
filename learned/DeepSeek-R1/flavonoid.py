"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:72010 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids have a 1-benzopyran (chromen-4-one) core with an aryl substituent at position 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define chromen-4-one core pattern (benzopyran-4-one)
    core_pattern = Chem.MolFromSmarts('c12c(ccc1)oc(=O)cc2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No chromen-4-one core found"
    
    # Find all core matches and check for aryl substituent
    matches = mol.GetSubstructMatches(core_pattern)
    for match in matches:
        try:
            # Get oxygen atom in the pyrone ring (index 3 in SMARTS pattern)
            oxygen_idx = match[3]
            oxygen = mol.GetAtomWithIdx(oxygen_idx)
            
            # Find adjacent carbon in pyrone ring (C2 position)
            pyrone_carbons = [n for n in oxygen.GetNeighbors() if n.GetSymbol() == 'C' and n.IsInRing()]
            if not pyrone_carbons:
                continue
                
            # Check substituents on each pyrone carbon adjacent to oxygen
            for carbon in pyrone_carbons:
                for bond in carbon.GetBonds():
                    neighbor = bond.GetOtherAtom(carbon)
                    # Look for aromatic substituents (aryl group)
                    if neighbor.GetIsAromatic():
                        # Verify substituent is part of an aromatic ring system
                        rings = mol.GetRingInfo().AtomRings()
                        for ring in rings:
                            if neighbor.GetIdx() in ring and len(ring) >= 6:
                                return True, "Chromen-4-one core with aryl substituent at position 2"
        except IndexError:
            continue
    
    return False, "No aryl substituent at position 2 of chromen-4-one core"
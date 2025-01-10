"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the isoflavonoid core (3-phenylchromen-4-one)
    # The SMARTS pattern includes atom mapping to identify position 3 and the aryl substituent
    isoflavonoid_smarts = """
    [O]=C1/C=C(\C=C2/COC=C12)[c:3][cH][cH][cH][cH][cH]
    """
    # Clean up the SMARTS pattern
    isoflavonoid_smarts = isoflavonoid_smarts.replace('\n', '').replace(' ', '')
    pattern = Chem.MolFromSmarts(isoflavonoid_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern for isoflavonoid core"
    
    # Find substructure matches of the isoflavonoid core
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No isoflavonoid core found"
    
    # Check each match for the aryl substituent at position 3
    for match in matches:
        # The index of the atom mapped as [c:3] corresponds to the carbon at position 3
        # Verify that this carbon is connected to an aryl group
        pos3_idx = match[pattern.GetSubstructMatch(pattern)[2]]  # Position of [c:3]
        pos3_atom = mol.GetAtomWithIdx(pos3_idx)
        
        # Get neighbors of the carbon at position 3 not part of the core
        core_atoms = set(match)
        neighbors = [nbr for nbr in pos3_atom.GetNeighbors() if nbr.GetIdx() not in core_atoms]
        
        # Check if any neighbor leads to an aryl group
        aryl_found = False
        for nbr in neighbors:
            # Use Chem.FindAtomEnvironmentOfRadiusN to find if an aromatic ring is connected
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=3, atomId=nbr.GetIdx())
            amap = {}
            submol = Chem.PathToSubmol(mol, env, atomMap=amap)
            aro_rings = submol.GetRingInfo().AtomRings()
            for ring in aro_rings:
                if all(submol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    aryl_found = True
                    break
            if aryl_found:
                break
        if aryl_found:
            return True, "Contains isoflavonoid core with aryl substituent at position 3"
        else:
            continue  # Check next match if available
    
    return False, "Isoflavonoid core found, but no aryl substituent at position 3"

# Example usage:
# smiles = 'Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O'  # Daidzein, an isoflavonoid
# result, reason = is_isoflavonoid(smiles)
# print(result, reason)
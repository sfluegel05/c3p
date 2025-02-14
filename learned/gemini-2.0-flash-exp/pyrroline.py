"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:35419 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a 5-membered ring with one nitrogen and one double bond

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a 5-membered ring with one nitrogen and 1 double bond, but without being part of a fused system.
    # also make sure that the N atom has only two bonds
    pyrroline_pattern = Chem.MolFromSmarts("[NX2;H0][CX3]=[CX2][CX2][CX2]")
    
    #Check for presence of the pattern
    if not mol.HasSubstructMatch(pyrroline_pattern):
      return False, "No pyrroline core structure found"
    
    
    #Ensure that it's not a pyrrole
    pyrrole_pattern = Chem.MolFromSmarts("[nX1][cX3]=[cX2][cX2]=[cX2]")
    if mol.HasSubstructMatch(pyrrole_pattern):
        return False, "It's a pyrrole not a pyrroline"

    # Check that there is no other double bonds in the 5 member ring
    ring_pattern = Chem.MolFromSmarts("[NX2,nX1;H0]1[CX3,cX2]~[CX2,cX2]~[CX2,cX2]~[CX2,cX2]~1")
    ring = mol.GetSubstructMatches(ring_pattern)
    
    if ring:
        for match in ring:
            submol = Chem.PathToSubmol(mol, match)
            db_count = 0
            for b in submol.GetBonds():
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    db_count +=1
            if db_count > 1:
                return False, "More than one double bond in the 5 membered ring"
            
    
    return True, "Pyrroline structure detected"
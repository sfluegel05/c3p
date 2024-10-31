from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylhydrazines(smiles: str):
    """
    Determines if a molecule is a phenylhydrazine (hydrazine with phenyl substituent).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phenylhydrazine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Look for phenylhydrazine substructure
    # Match both N-N single bonds and N=N double bonds connected to phenyl
    phenylhydrazine_pattern = Chem.MolFromSmarts('[c]!@[N]-[N,n]')
    phenylhydrazine_pattern2 = Chem.MolFromSmarts('[c]!@[N]=[N,n]')
    phenylhydrazine_pattern3 = Chem.MolFromSmarts('[c]!@[N](-[*])-[N,n]')
    
    matches = []
    matches.extend(mol.GetSubstructMatches(phenylhydrazine_pattern))
    matches.extend(mol.GetSubstructMatches(phenylhydrazine_pattern2))
    matches.extend(mol.GetSubstructMatches(phenylhydrazine_pattern3))
    
    if not matches:
        return False, "No phenylhydrazine group found"

    # Verify the phenyl group is a proper 6-membered aromatic ring
    for match in matches:
        phenyl_atom = mol.GetAtomWithIdx(match[0])
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if phenyl_atom.GetIdx() in ring and len(ring) == 6:
                atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if all(atom.GetIsAromatic() and atom.GetSymbol() == 'C' for atom in atoms):
                    return True, "Contains phenylhydrazine group"
                    
    return False, "No proper phenyl ring connected to hydrazine"
# Pr=None
# Recall=None
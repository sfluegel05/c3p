"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:24768 enone
Definition: An alpha,beta-unsaturated ketone of general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) ≠ H) 
in which the C=O function is conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha,beta-unsaturated ketone substructure
    enone_pattern = Chem.MolFromSmarts("[CD3]=[CD2][C@]=O")
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    
    if not enone_matches:
        return False, "No alpha,beta-unsaturated ketone substructure found"
    
    # Check for conjugation and non-hydrogen substituent at carbonyl carbon
    for match in enone_matches:
        c_alpha, c_beta, c_carbonyl = match
        
        # Check for non-hydrogen substituent at carbonyl carbon
        carbonyl_atom = mol.GetAtomWithIdx(c_carbonyl)
        if all(nbr.GetAtomicNum() == 1 for nbr in carbonyl_atom.GetNeighbors()):
            continue
        
        # Check if alpha,beta-unsaturated ketone is part of an aromatic system
        aromatic_rings = mol.GetAromaticRings()
        for ring in aromatic_rings:
            if c_alpha in ring and c_beta in ring and c_carbonyl in ring:
                return True, "Molecule contains an alpha,beta-unsaturated ketone with conjugation"
    
    return False, "Alpha,beta-unsaturated ketone not conjugated or no non-hydrogen substituent at carbonyl carbon"
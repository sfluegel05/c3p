from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_furoate_ester(smiles: str):
    """
    Determines if a molecule is a 2-furoate ester (carboxylic ester where the carboxylic acid component is 2-furoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-furoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for 2-furoate ester group:
    # [#6]-O-C(=O)-c1ccco1
    # Where:
    # [#6] is any carbon (the alcohol part)
    # -O-C(=O)- is the ester linkage
    # c1ccco1 is the 2-furan ring
    furoate_pattern = Chem.MolFromSmarts('[#6]-O-C(=O)-c1ccco1')
    
    if mol.HasSubstructMatch(furoate_pattern):
        # Get the matches to analyze the furan ring
        matches = mol.GetSubstructMatches(furoate_pattern)
        
        for match in matches:
            # The last 5 atoms in the match correspond to the furan ring
            furan_atoms = match[-5:]
            furan_ring_atoms = [mol.GetAtomWithIdx(i) for i in furan_atoms]
            
            # Verify it's a proper furan ring (aromatic and has oxygen)
            if all(atom.GetIsAromatic() for atom in furan_ring_atoms) and \
               any(atom.GetSymbol() == 'O' for atom in furan_ring_atoms):
                
                # Get the alcohol component (first atom in match)
                alcohol_carbon = mol.GetAtomWithIdx(match[0])
                
                return True, "Contains 2-furoate ester group"
                
    return False, "Does not contain 2-furoate ester group"
# Pr=1.0
# Recall=1.0
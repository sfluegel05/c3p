from rdkit import Chem
from rdkit.Chem import AllChem

def is_monohydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a monohydroxyflavanone (flavanone with exactly one hydroxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monohydroxyflavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavanone core structure
    flavanone_pattern = Chem.MolFromSmarts('[O;R1]-1-c2ccccc2-C(=O)C[C@H]-1c1ccccc1')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone structure"

    # Count number of hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)

    if num_hydroxy == 0:
        return False, "No hydroxy groups found"
    elif num_hydroxy > 1:
        return False, f"Found {num_hydroxy} hydroxy groups - more than one"
    else:
        # Get position of hydroxy group
        hydroxy_atom = mol.GetAtomWithIdx(hydroxy_matches[0][0])
        position = []
        
        # Try to determine position of OH group
        for match in mol.GetSubstructMatches(flavanone_pattern):
            core_atoms = set(match)
            for atom in mol.GetAtomWithIdx(hydroxy_matches[0][0]).GetNeighbors():
                if atom.GetIdx() in core_atoms:
                    position.append(atom.GetIdx())
                    
        if position:
            return True, f"Monohydroxyflavanone with OH at position {position[0]}"
        else:
            return True, "Monohydroxyflavanone"
# Pr=1.0
# Recall=0.72
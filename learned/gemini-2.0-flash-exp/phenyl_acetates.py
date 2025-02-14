"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate has an acetate group attached to a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phenyl ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"
    
    # Check for acetate group
    acetate_pattern = Chem.MolFromSmarts("C(=O)OC")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    if len(acetate_matches) == 0:
            return False, "No acetate group found"
    
    #Verify attachment of acetate to phenyl
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)

    found_attachment = False
    for acetate_match in acetate_matches:
        for phenyl_match in phenyl_matches:
           # Get the oxygen of acetate (second to last)
            oxygen_index = acetate_match[-2]
            acetate_oxygen_atom = mol.GetAtomWithIdx(oxygen_index)
            #Check the neighbors of the oxygen and if a neighbor is in the phenyl ring
            for neighbor in acetate_oxygen_atom.GetNeighbors():
                if neighbor.GetIdx() in phenyl_match:
                    found_attachment = True
                    break
            if found_attachment:
                break
        if not found_attachment:
            return False, "Acetate is not directly attached to the phenyl ring"
        else:
            found_attachment=False #reset to find another match
            
    
    
    return True, "Contains a phenyl ring with at least one acetate group directly attached"
"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with a single acyl, alkyl, or alk-1-enyl
    substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify glycerol backbone pattern: C-C-C with at least two hydroxyls
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Collect substituents attached to the backbone
    long_chain_smarts = "[R]C(=O)C([R])"  # Common long chains in lipids
    potential_substituents = mol.GetSubstructMatches(Chem.MolFromSmarts(long_chain_smarts))
    
    # Further filter to ensure connection to one glycerol center only
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    # We assume all substituents are connected to one glycerol carbon
    num_glycerol_substituents = 0
    
    for a, b, c in glycerol_matches:
        # Get the neighbors of the central glycerol carbon that are part of larger chains
        # counting substituents directly connected to glycerol carbons
        central_carbon_conn = mol.GetAtomWithIdx(a+1).GetNeighbors()
        substituent_atoms = [atom for atom in central_carbon_conn if atom.GetIdx() not in [a, b, c]]
        
        for sub_atom in substituent_atoms:
            for substituent in potential_substituents:
                if sub_atom.GetIdx() in substituent:
                    num_glycerol_substituents += 1

    if num_glycerol_substituents != 1:
        return False, f"Expected one connected substituent, found {num_glycerol_substituents}"
    
    return True, "Contains glycerol backbone with exactly one lipid substituent"
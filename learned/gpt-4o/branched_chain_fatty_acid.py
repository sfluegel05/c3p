"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A BCFA is characterized by the presence of a carboxylic acid group (or its ionized form),
    and one or more branching points typically involving smaller alkyl groups attached 
    to the main hydrocarbon chain, avoiding ring-contained branches except cyclopropyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: A tuple containing a boolean value indicating if it's a BCFA, and the reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for carboxylic acid group (including ionized form)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify branching points: carbon with 3 or more non-ring single bonds
    branch_found = False
    for atom in mol.GetAtoms():
        # Check if the atom is carbon
        if atom.GetAtomicNum() != 6:
            continue
        # Check connections and ensure it is not part of a ring
        if atom.GetDegree() > 2 and not atom.IsInRing():
            branch_found = True
            break
            
    if not branch_found:
        return False, "No suitable branching pattern found"

    return True, "Contains a branched hydrocarbon chain with carboxylic acid group"

# Test the function with some example SMILES strings
smiles_examples = [
    "OC(=O)\\C=C\\C(C)(C)C",  # 4,4-dimethyl-2E-pentenoic acid
    "OC(=O)CCC(C)(C)C",      # 4,4-dimethyl valeric acid
    "CC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylnonadecanoic acid
    "CC(C)(C)C(O)=O",        # pivalic acid
    "CCCCCCCCCCCCCCCCCC(=O)O",  # Straight chain - not branched
]

# Output results for testing examples
for sm in smiles_examples:
    result, reason = is_branched_chain_fatty_acid(sm)
    print(f"SMILES: {sm} -> Result: {result}, Reason: {reason}")
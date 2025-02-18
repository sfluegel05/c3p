"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: Phosphatidyl-L-serine, a class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine is defined as an aminophospholipid in which a phosphatidyl group is esterified
    to the hydroxy group of serine. In our heuristic, the molecule must contain:
      - at least one phosphorus atom,
      - a phosphoserine head-substructure (using a SMARTS for a phosphate group attached to OCC(N)C(=O)O),
      - and at least 2 long alkyl chain fragments (defined here as a chain of 6 consecutive carbons) 
        as surrogates for fatty acyl chains.
        
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as phosphatidyl-L-serine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of phosphorus atoms.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain any phosphorus atoms"
    
    # Define a SMARTS pattern for the phosphoserine head group.
    # We look for a phosphate (P(=O)(O)(...)) where one substituent is OCC(N)C(=O)O.
    phosphoserine_smarts = "P(=O)(O)(OCC(N)C(=O)O)"
    phosphoserine_pattern = Chem.MolFromSmarts(phosphoserine_smarts)
    if phosphoserine_pattern is None:
        return False, "Error creating phosphoserine SMARTS pattern"
    
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Phosphoserine head group not found"
    
    # Define a SMARTS pattern for a long aliphatic chain fragment (at least 6 carbons).
    # This is a heuristic for a fatty acyl chain.
    alkyl_chain_smarts = "CCCCCC"
    alkyl_chain_pattern = Chem.MolFromSmarts(alkyl_chain_smarts)
    if alkyl_chain_pattern is None:
        return False, "Error creating alkyl chain SMARTS pattern"
    
    chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(chain_matches) < 2:
        return False, f"Found only {len(chain_matches)} long alkyl chain fragments; expected at least 2 for diacyl phospholipid"
    
    # Optional additional checks: molecular weight (PS tend to be heavy lipid molecules).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low to be a phosphatidyl-L-serine lipid"
    
    return True, "Molecule contains phosphoserine head group and sufficient acyl chains to be phosphatidyl-L-serine"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES examples:
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"
    result, reason = is_phosphatidyl_L_serine(test_smiles)
    print(result, "->", reason)
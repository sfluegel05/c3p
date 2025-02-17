"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: Phosphatidylglycerol 
Definition: A glycerophosphoglycerol that is glycerol in which the hydrogen of one of 
the primary hydroxy groups has been replaced by a phosphatidyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    A phosphatidylglycerol has a single phosphorus atom that is bound in a phosphoglycerol headgroup 
    (i.e. P(=O)(O)(OCC(O)CO)) and contains two acyl chains (detected as ester groups "OC(=O)") attached 
    to a glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that there is exactly one phosphorus atom (atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 1:
        return False, f"Expected exactly one phosphorus atom, found {len(p_atoms)}"
    
    # Define a SMARTS for the phosphoglycerol headgroup.
    # This pattern looks for a phosphorus double-bonded to an oxygen and having one substituent that 
    # is an oxygen bound to a 2-carbon chain bearing two hydroxyl groups: OCC(O)CO.
    pg_head_smarts = "P(=O)(O)(OCC(O)CO)"
    pg_head = Chem.MolFromSmarts(pg_head_smarts)
    if not mol.HasSubstructMatch(pg_head):
        return False, "Phosphoglycerol headgroup pattern not found"
    
    # Define a SMARTS for acyl ester groups found in diacylglycerol parts ("OC(=O)")
    ester_smarts = "OC(=O)"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # For phosphatidylglycerol we expect two acyl chains (two ester bonds) attached to the glycerol backbone.
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} acyl ester(s); expected 2 acyl chains"
    
    # Optionally, one can check with RDKit's rotatable bonds or molecular weight if desired.
    # Here we simply rely on the presence of the required substructures.
    
    return True, "Molecule contains a phosphoglycerol headgroup and two acyl ester chains consistent with phosphatidylglycerol"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: PG(8:0/8:0)
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OC[C@@H](O)CO)(O)=O"
    result, reason = is_phosphatidylglycerol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)